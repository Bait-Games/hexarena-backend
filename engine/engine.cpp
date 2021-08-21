#include "engine.h"
#include "utils.h"

#include <cmath>
#include <random>
#include <iostream>
#include <thread>

/**
 * Starts the engine in a new thread.
 */
Engine::Engine() {
    this->tick_num = 0;
    this->next_user_id = 0;
    this->running = true;
    this->runner = std::thread(&Engine::run, this);
}

/**
 * Adds a new user at the given position.
 * @param x The x coordinate of the user's position.
 * @param y The y coordinate of the user's position.
 * @return The id of the user
 */
int Engine::add_user(double x, double y) {
    int id = next_user_id++;
    return this->users.insert({id++, {id, x, y, [&](int amt, double x, double y){
        this->drop(amt, x, y);
    }}}).first->first;
}

/**
 * Moves the given user towards the direction at the given angle.
 * @param id The id of the user to move.
 * @param angle The angle of the vector pointing towards the desired movement direction.
 */
void Engine::move(int id, double angle) {
    std::lock_guard<std::mutex> guard(this->mutex);
    User &user = this->users.at(id);
    if (user.moved) return;
    user.moved = true;
    user.momentum_x += std::cos(angle) * MOMENTUM_GAIN_PER_TICK;
    user.momentum_y += std::sin(angle) * MOMENTUM_GAIN_PER_TICK;
}

/**
 * Rotates the given user to face towards the given angle
 * @param id The id of the user to rotate.
 * @param angle The angle of the vector pointing towards the desired facing.
 */
void Engine::rotate(int id, double angle) {
    std::lock_guard<std::mutex> guard(this->mutex);
    User &user = this->users.at(id);
    if (user.rotated) return;
    user.rotated = true;

    double size = user.get_size();
    for (int iter = 0; user.rotation != angle && iter < MAX_ROTATION_ITERATIONS; iter++) {
        double direction = angle > user.rotation && angle < user.rotation + M_PI ? 1 : -1;

        double delta = MAX_SAFE_MOVEMENT / size * direction;
        double new_angle = user.rotation + delta;

        if ((user.rotation < angle && angle < new_angle) || (user.rotation > angle && angle > new_angle)) {
            new_angle = angle;
        }

        double old_angle = user.rotation;
        user.rotation = new_angle;
        bool blocked = false;
        // user2 is a copy, do not attempt to modify it or pass it anywhere as a reference
        for (auto& [id2, user2] : this->users) {
            if (id == id2) continue;
            blocked = blocked || user.collide_with_user(this->users.at(id2));
        }

        if (blocked) {
            user.rotation = old_angle;
            break;
        }
    }
}

/**
 * Creates a new user.
 * @return The id of the new user.
 */
int Engine::create() {
    std::lock_guard<std::mutex> guard(this->mutex);
    std::uniform_real_distribution<double> uniform_width(0, WORLD_WIDTH);
    std::uniform_real_distribution<double> uniform_height(0, WORLD_HEIGHT);
    double x, y;

    bool has_space;
    int iter = 0;
    do {
        x = uniform_width(gen);
        y = uniform_height(gen);
        has_space = true;
        for (auto [id, user] : users) {
            if (distance(x, y, user.x, user.y) < user.get_size() + CELL_INNER_RADIUS * 2) {
                has_space = false;
            }
        }
    } while (!has_space && iter++ < 20);

    if (!has_space) {
        std::cerr << "no space for new user" << std::endl;
        return -1;
    }

    return this->add_user(x, y);
}

/**
 * Removes (kills) the given user.
 * @param id The id of the user to remove.
 */
void Engine::remove(int id) {
    std::lock_guard<std::mutex> guard(this->mutex);
    User &user = this->users.at(id);
    user.shrink(0);
    user.destroyed = true;
    if (user.resources > 0) {
        this->resources.push_back({user.x, user.y, user.resources / 2});
    }
    for (auto observer : user.observers) // TODO: discuss how to notify observers of destruction;
    this->users.erase(id);
}

/**
 * Attaches a new bodypart to the given user.
 * @param id The id of the user to modify.
 * @param type The type of the bodypart to add.
 * @param part The id of the cell to which the new part should be attached.
 * @param face The face of the cell at which the new bodypart should be attached.
 * @return 0 on success, error code on failure:
 *   - -1 - asked to attach a component to a component that is not a cell
 *   - -2 - asked to attach a component on a face that is not free
 *   - -3 - user with given id does not exist.
 *   - -4 - the given face of the given part is free but another cell occupies the space in front of it.
 *          This should never happen, but is controlled as a sanity check.
 *   - -5 - asked to attach a component to a component that does not exist
 *   - -6 - trying to build a cell on a space that is already occupied by a (small) component
 *   - -7 - user does not have sufficient resources to buy requested bodypart.
 *   - -8 - trying to build a component of type NONE
 */
int Engine::attach(int id, BODYPART_TYPE type, int part, int face) {
    std::lock_guard<std::mutex> guard(this->mutex);
    if (!this->users.contains(id)) return -3;
    User &user = this->users.at(id);
    if (user.resources < BODYPART_COST.at(type)) return -7;

    int ret = user.grow(part, face, type);
    if (ret == 0) user.resources -= BODYPART_COST.at(type);
    return ret;
}

/**
 * Removes the given bodypart from the given user.
 * @param id The id of the user to modify.
 * @param part The id of the bodypart to remove.
 * @return true iff this caused the user's death, i.e. the removed bodypart was the root.
 */
int Engine::detach(int id, int part) {
    std::lock_guard<std::mutex> guard(this->mutex);
    if (!this->users.contains(id)) return -3;
    User &user = this->users.at(id);
    if (user.next_component > part || user.components[part].type == NONE) return -2;

    return user.shrink(part);
}

/**
 * Gets information about a given user.
 * @param id The id of the user to retrieve.
 * @return The User object of the given user.
 */
const User& Engine::user(int id) {
    std::lock_guard<std::mutex> guard(this->mutex);
    return this->users.at(id);
}
const User& Engine::info(int id) {
    return this->user(id);
}

/**
 * Adds the given observer to the world. It will be called every tick and passed the current state of the world.
 * @param cb The callback to be notified every tick.
 */
void Engine::observe(void (*cb)(const std::vector<Resource> &, const std::map<int, User> &)) {
    this->observers.push_back(cb);
}
void Engine::register_global(void (*cb)(const std::vector<Resource> &, const std::map<int, User> &)) {
    this->observe(cb);
}

/**
 * Adds the given observer to the given user. It will be called evert tick and passed the current state of the user.
 * @param user The user to observe.
 * @param cb The callback to be notified every tick.
 */
void Engine::observe(int user, void (*cb)(const User &)) {
    this->users.at(user).observers.push_back(cb);
}

/**
 * Creates a resource packet of the given size at the given coordinates
 * @param amt The amount of resources to drop.
 * @param x The x coordinate of the resource drop's position.
 * @param y The y coordinate of the resource drop's position.
 */
void Engine::drop(int amt, double x, double y) {
    this->resources.push_back({x, y, amt});
}

/**
 * Processes everything that goes on in a single tick
 */
void Engine::tick() {
    std::lock_guard<std::mutex> guard(this->mutex);
    this->tick_num++;

    if (this->resources.size() < RESOURCE_DENSITY * WORLD_WIDTH * WORLD_HEIGHT / (1024*1024)) {
        std::uniform_real_distribution<double> uniform_width(0, WORLD_WIDTH);
        std::uniform_real_distribution<double> uniform_height(0, WORLD_HEIGHT);
        this->resources.push_back({uniform_width(gen), uniform_height(gen), NATURAL_RESOURCE_AMOUNT});
    }

    for (auto [id, _] : this->users) {
        User &user = this->users.at(id);
        user.tick_reset();
        user.tick_parts();

        double dx = user.momentum_x;
        double dy = user.momentum_y;
        while (dx > 0 || dy > 0) {
            double step_x = dx, step_y = dy;

            double len = distance(0, 0, dx, dy);
            if (len > MAX_SAFE_MOVEMENT) {
                double frac = MAX_SAFE_MOVEMENT / len;
                step_x = frac * dx;
                step_y = frac * dy;
            }
            dx -= step_x;
            dy -= step_y;

            user.x += step_x;
            user.y += step_y;
            user.momentum_x *= MOMENTUM_MULTIPLIER_PER_TICK;
            user.momentum_y *= MOMENTUM_MULTIPLIER_PER_TICK;

            bool blocked = false;
            for (auto [id2, _2] : this->users) {
                if (id == id2) continue;
                if (user.collide_with_user(this->users.at(id2))) blocked = true;
            }

            for (Resource &res : this->resources) {
                if (user.collide_with_resource(res)) blocked = true;
            }

            if (blocked) {
                user.x -= step_x;
                user.y -= step_y;
                break;
            }
        }

        if (user.y > WORLD_HEIGHT) user.y = WORLD_HEIGHT;
        if (user.y < 0) user.y = 0;
        if (user.x > WORLD_WIDTH) user.x = WORLD_WIDTH;
        if (user.x < 0) user.x = 0;
    }

    std::erase_if(this->resources, [](const Resource& res){
        return res.amt <= 0;
    });

    std::erase_if(this->users, [](const User &user){
        // TODO: observers need to be notified of death.
        return user.destroyed;
    });

    for (auto observer : this->observers) {
        observer(this->resources, this->users);
    }
    for (auto [id, user] : this->users) {
        for (auto observer : user.observers) {
            observer(user);
        }
    }
}
/**
 * Loops until running is false, calling tick every 1/TICK_RATE seconds.
 */
void Engine::run() {
    std::chrono::time_point<std::chrono::steady_clock> last_tick = std::chrono::steady_clock::now();
    while (this->running) {
        uint time_between_ticks = 1000 / TICK_RATE;
        std::this_thread::sleep_until(last_tick + std::chrono::milliseconds(time_between_ticks));
        last_tick = std::chrono::steady_clock::now();
        tick();
    }
}

/**
 * Halts the engine.
 */
Engine::~Engine() {
    this->running = false;
    this->runner.join();
}
