#include "engine.h"
#include "utils.h"

#include <cmath>
#include <random>
#include <iostream>
#include <thread>

Engine::Engine() {
    this->tick_num = 0;
    this->next_user_id = 0;
    this->running = true;
    this->runner = std::thread(&Engine::run, this);
}

int Engine::add_user(double x, double y) {
    int id = next_user_id++;
    return this->users.insert({id++, {id, x, y, [&](int amt, double x, double y){
        this->drop(amt, x, y);
    }}}).first->first;
}

void Engine::move(int id, double angle) {
    std::lock_guard<std::mutex> guard(this->mutex);
    User &user = this->users.at(id);
    if (user.moved) return;
    user.moved = true;
    user.momentum_x += std::cos(angle) * MOMENTUM_GAIN_PER_TICK;
    user.momentum_y += std::sin(angle) * MOMENTUM_GAIN_PER_TICK;
}

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

int Engine::attach(int id, BODYPART_TYPE type, int part, int face) {
    std::lock_guard<std::mutex> guard(this->mutex);
    if (!this->users.contains(id)) return -3;
    User &user = this->users.at(id);
    if (user.resources < BODYPART_COST.at(type)) return -7;

    int ret = user.grow(part, face, type);
    if (ret == 0) user.resources -= BODYPART_COST.at(type);
    return ret;
}

int Engine::detach(int id, int part) {
    std::lock_guard<std::mutex> guard(this->mutex);
    if (!this->users.contains(id)) return -3;
    User &user = this->users.at(id);
    if (user.next_component > part || user.components[part].type == NONE) return -2;

    return user.shrink(part);
}

const User& Engine::user(int id) {
    std::lock_guard<std::mutex> guard(this->mutex);
    return this->users.at(id);
}
const User& Engine::info(int id) {
    return this->user(id);
}

void Engine::observe(void (*cb)(const std::vector<Resource> &, const std::map<int, User> &)) {
    this->observers.push_back(cb);
}
void Engine::register_global(void (*cb)(const std::vector<Resource> &, const std::map<int, User> &)) {
    this->observe(cb);
}

void Engine::observe(int user, void (*cb)(const User &)) {
    this->users.at(user).observers.push_back(cb);
}

void Engine::drop(int amt, double x, double y) {
    this->resources.push_back({x, y, amt});
}

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
void Engine::run() {
    std::chrono::time_point<std::chrono::steady_clock> last_tick = std::chrono::steady_clock::now();
    while (this->running) {
        uint time_between_ticks = 1000 / TICK_RATE;
        std::this_thread::sleep_until(last_tick + std::chrono::milliseconds(time_between_ticks));
        last_tick = std::chrono::steady_clock::now();
        tick();
    }
}

Engine::~Engine() {
    this->running = false;
    this->runner.join();
}
