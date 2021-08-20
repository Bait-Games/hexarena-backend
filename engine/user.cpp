#include <iostream>
#include <cmath>
#include <utility>

#include "user.h"

User::User(int id, double x, double y, std::function<void(int, double, double)> drop) {
    this->id = id;
    this->x = x;
    this->y = y;
    this->momentum_x = 0;
    this->momentum_y = 0;
    this->rotation = 0;
    this->resources = 0;
    this->kills = 0;
    this->components = new Component[16];
    this->components_size = 16;
    for (Component* it = this->components; it < this->components + this->components_size; it++) {
        it->type = NONE;
    }
    this->next_component = 0;
    this->destroyed = false;
    this->moved = false;
    this->rotated = false;
    this->size = CELL_INNER_RADIUS * 2;
    this->drop = std::move(drop);
}

void User::tick_reset() {
    moved = false;
    rotated = false;
}

void User::tick_parts() const {
    for (Component* it = this->components; it < this->components + this->components_size; it++) {
        switch(it->type) {
            case CELL:
                it->health += REGEN_RATE;
                if (it->health > MAX_HEALTH) it->health = MAX_HEALTH;
                break;
            case BOUNCE:
                it->inflated += INFLATE_RATE;
                if (it->inflated >= MAX_INFLATE) it->inflated = MAX_INFLATE;
            case SPIKE: // NOLINT(bugprone-branch-clone)
            case SHIELD:
                break;
            case NONE:
                // skip - this cell of the array is empty
                break;
        }
    }
}


double distance(User &a, User &b) {
    return distance(a.x, a.y, b.x, b.y);
}
std::pair<double, double> rel_pos(double x, double y, int up, int fwd, int bwd) {
#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnusedValue"
    y -= fwd * CELL_INNER_RADIUS * 2;
    bwd += fwd;
    fwd -= fwd;

    // now up and bwd are opposites by the up+fwd+bwd=0 equality
    y -= CELL_INNER_RADIUS * up;
    x += CELL_INNER_RADIUS + 20 * up;

    return {x, y};
#pragma clang diagnostic pop
}

/**
 * Attaches a new part to the user. Specifying an invalid face (< 0 || > 6) is undefined behavior
 * @param part The id of the component at which the new component is attached
 * @param face The face on the component at which the new component is attached
 * @param type The type of the component to attach
 * @return 0 on success, error code on failure:
 *   - -1 - asked to attach a component to a component that is not a cell
 *   - -2 - asked to attach a component on a face that is not free
 *   - -3 - reserved
 *   - -4 - the given face of the given part is free but another cell occupies the space in front of it.
 *          This should never happen, but is controlled as a sanity check.
 *   - -5 - asked to attach a component to a component that does not exist
 *   - -6 - trying to build a cell on a space that is already occupied by a (small) component
 *   - -7 - reserved
 *   - -8 - trying to build a component of type NONE
 */
int User::grow(int part, int face, BODYPART_TYPE type)  {
    if (type == NONE) return -8;
    if (part < 0 || part >= next_component || this->components[part].type == NONE) return -5;
    if (this->components[part].type != CELL) return -1;
    if (this->components[part].faces[face] != -1) return -2;

    int newcoords[3] = {
            std::get<0>(this->components[part].coords),
            std::get<1>(this->components[part].coords),
            std::get<2>(this->components[part].coords)
    };
    if (face == 2 || face == 3) {
        newcoords[0] += 1;
    }
    if (face == 5 || face == 0) {
        newcoords[0] -= 1;
    }
    if (face == 0 || face == 1) {
        newcoords[1] += 1;
    }
    if (face == 3 || face == 4) {
        newcoords[1] -= 1;
    }
    if (face == 4 || face == 5) {
        newcoords[2] += 1;
    }
    if (face == 1 || face == 2) {
        newcoords[2] -= 1;
    }

    // if everything else works this should be impossible. Once we're convinced it is so we can remove this check
    auto iter = std::find_if(this->components, (this->components + this->next_component),
                 [&](Component c){
        return c.type == CELL
        && std::get<0>(c.coords) == newcoords[0]
        && std::get<1>(c.coords) == newcoords[1]
        && std::get<2>(c.coords) == newcoords[2];
    });
    if (iter < this->components + this->next_component) {
        std::cerr << "inconsistent components: face free but space occupied" << std::endl;
        return -4;
    }

    // check that no other bodypart (particularly spikes/shields/bouncers) occupies the space we're trying to fill
    if (type == CELL) {
        iter = std::find_if(this->components, this->components + this->next_component,
                  [&](Component c){
                      if (c.type != CELL) return false;
                      // returns true if there's a cell adjacent to our target position with the face pointing to us
                      // occupied. That would mean there's a component (e.g. a spike) where we're trying to build.
                      return std::get<0>(c.coords) == newcoords[0] - 1
                      && std::get<1>(c.coords) == newcoords[1] + 1
                      && std::get<2>(c.coords) == newcoords[2]
                      && c.faces[5] != -1
                      ||
                      std::get<0>(c.coords) == newcoords[0]
                      && std::get<1>(c.coords) == newcoords[1] + 1
                      && std::get<2>(c.coords) == newcoords[2] - 1
                      && c.faces[4] != -1
                      ||
                      std::get<0>(c.coords) == newcoords[0] + 1
                      && std::get<1>(c.coords) == newcoords[1]
                      && std::get<2>(c.coords) == newcoords[2] - 1
                      && c.faces[3] != -1
                      ||
                      std::get<0>(c.coords) == newcoords[0] + 1
                      && std::get<1>(c.coords) == newcoords[1] - 1
                      && std::get<2>(c.coords) == newcoords[2]
                      && c.faces[2] != -1
                      ||
                      std::get<0>(c.coords) == newcoords[0]
                      && std::get<1>(c.coords) == newcoords[1] - 1
                      && std::get<2>(c.coords) == newcoords[2] + 1
                      && c.faces[1] != -1
                      ||
                      std::get<0>(c.coords) == newcoords[0] - 1
                      && std::get<1>(c.coords) == newcoords[1]
                      && std::get<2>(c.coords) == newcoords[2] + 1
                      && c.faces[0] != -1;
        });
        if (iter < this->components + this->next_component) return -6;
    }

    if (next_component == components_size) {
        components_size *= 2;
        auto* tmp = new Component[components_size];
        for (int i = 0; i < this->next_component; i++) {
            tmp[i] = this->components[i];
        }
        delete this->components;
        this->components = tmp;
    }

    int i = this->next_component++;
    this->components[i].type = type;

    switch (type) {
        case BOUNCE:
            this->components[i].inflated = MAX_INFLATE;
        case SPIKE:
        case SHIELD:
            this->components[i].body = part;
            break;
        case CELL:
            this->components[i].health = MAX_HEALTH;
            this->components[i].dead = false;
            for (int & _face : this->components[i].faces) _face = -1;
            this->components[i].coords = {newcoords[0], newcoords[1], newcoords[2]};

            // set the faces for all adjacent cells
            for (int j = 0; j < this->next_component - 1; j++) {
                Component* c = this->components + j;
                if (c->type != CELL) continue;
                if (c->coords == (std::tuple<int, int, int>){newcoords[0]-1, newcoords[1]+1, newcoords[2]}) {
                    this->components[i].faces[2] = j;
                    c->faces[5] = i;
                } else if (c->coords == (std::tuple<int, int, int>){newcoords[0], newcoords[1]+1, newcoords[2]-1}) {
                    this->components[i].faces[1] = j;
                    c->faces[4] = i;
                } else if (c->coords == (std::tuple<int, int, int>){newcoords[0]+1, newcoords[1], newcoords[2]-1}) {
                    this->components[i].faces[0] = j;
                    c->faces[3] = i;
                } else if (c->coords == (std::tuple<int, int, int>){newcoords[0]+1, newcoords[1]-1, newcoords[2]}) {
                    this->components[i].faces[5] = j;
                    c->faces[2] = i;
                } else if (c->coords == (std::tuple<int, int, int>){newcoords[0], newcoords[1]-1, newcoords[2]+1}) {
                    this->components[i].faces[4] = j;
                    c->faces[1] = i;
                } else if (c->coords == (std::tuple<int, int, int>){newcoords[0]-1, newcoords[1], newcoords[2]+1}) {
                    this->components[i].faces[3] = j;
                    c->faces[0] = i;
                }
            }
            break;
        case NONE:
            // cannot occur, we made sure type is not none earlier
            break;
    }
    this->components[part].faces[face] = this->next_component - 1;
    this->compute_size();
    return 0;
}

/**
 * Damages the given component by the given amount.
 * @param part the component to damage. Must be of type cell.
 * @param amt the amount of damage taken.
 * @return true iff the component was destroyed
 */
bool User::damage(int part, int amt) {
    Component* component = this->components + part;

    if (component->type == SHIELD) return false;
    if (component->type == SPIKE || component->type == BOUNCE) component = this->components + component->body;
    // now component has to point to a cell
    component->health -= amt;
    if (component->health <= 0) {
        // cast is safe - we will never have more than MAX_INT components
        return this->shrink(static_cast<int>(std::distance(this->components, components)));
    }
    return false;
}

/**
 * Destroys the given component of this user.
 * @param part The id (index in array) of the component to destroy
 * @return true iff this shrinkage caused user death
 */
bool User::shrink(int part) {
    if (part == 0) {
        if (this->components[0].dead) return false;
        this->components[0].dead = true;
        if (!this->destroyed) {
            this->destroyed = true;
            if (this->resources > 0) {
                this->drop(this->resources/2, this->x, this->y);
            }
        }
        for (int neighbor : this->components[0].faces) {
            if (neighbor == -1) continue;
            this->shrink(neighbor);
        }
        return true;
    }
    if (this->components[part].type == NONE) return false;
    if (this->components[part].type == CELL) {
        auto [x, y, _] = this->get_hitbox(part);
        this->drop(BODYPART_COST.at(CELL)/2, x, y);
    }
    this->components->type = NONE;
    for (int i = 0; i < this->next_component; i++) {
        if (this->components[i].type != CELL) continue;
        for (int &face : this->components[i].faces) {
            if (face == part) face = -1;
        }
    }
    std::set<int> connected;
    this->find_connected(0, connected);
    for (int i = 0; i < this->next_component; i++) {
        if (connected.contains(i)) continue;
        if (this->components[i].type == CELL) {
            auto [x, y, _] = this->get_hitbox(i);
            this->drop(BODYPART_COST.at(CELL)/2, x, y);
        }
        this->components[i].type = NONE;
    }
    this->compute_size();
    return false;
}

void User::find_connected(int node, std::set<int> &set) {
    if (set.contains(node)) return;
    set.insert(node);
    if (this->components[node].type != CELL) return;
    for (int face : this->components[node].faces) {
        if (face == -1) continue;
        this->find_connected(face, set);
    }
}

/**
 * Given the id (index in array) of a component of this user, returns its current hitbox.
 * @param component The id of the component. Index into this->components. Must be CELL or SPIKE
 * @return a tuple (x, y, size) where x and y are the coordinates of the center of the hitbox and size is its radius.
 *         All hitboxes are circular
 */
std::tuple<double, double, int> User::get_hitbox(int component) const {
    double pos_x, pos_y;
    int size;
    if (this->components[component].type == CELL) {
        auto [pos_x_tmp, pos_y_tmp] = rel_pos(this->x, this->y,
                                              std::get<0>(this->components[component].coords),
                                              std::get<1>(this->components[component].coords),
                                              std::get<2>(this->components[component].coords));
        pos_x = pos_x_tmp;
        pos_y = pos_y_tmp;
        size = CELL_INNER_RADIUS;
    } else {
        int parent = this->components[component].body;
        auto [parent_x, parent_y] = rel_pos(this->x, this->y,
                                            std::get<0>(this->components[parent].coords),
                                            std::get<1>(this->components[parent].coords),
                                            std::get<2>(this->components[parent].coords));
        int face = -1;
        for (int i = 0; i < 6; i++) {
            if (this->components[parent].faces[i] == component) {
                face = i;
                break;
            }
        }
        int mult_x, mult_y;
        switch (face) {
            case 0:
                mult_x = -24;
                mult_y = -14;
                break;
            case 1:
                mult_x = 0;
                mult_y = -28;
                break;
            case 2:
                mult_x = 24;
                mult_y = -14;
                break;
            case 3:
                mult_x = 24;
                mult_y = 14;
                break;
            case 4:
                mult_x = 0;
                mult_y = 28;
                break;
            case 5:
                mult_x = -24;
                mult_y = 14;
                break;
            default:
                std::cerr << "Malformed body: parent (" << parent << ") not adjacent to child (" << component << ")" << std::endl;
                throw std::invalid_argument("Malformed body: parent not adjacent to child");
        }

        switch (this->components[component].type) {
            // We do not use hitboxes for shield and bounce. The only things that need hitbox detection are spikes and cells
//            case BOUNCE:
//                pos_x = parent_x + (10./14) * mult_x;
//                pos_y = parent_y + (10./14) * mult_y;
//                size = 18;
//                break;
//            case SHIELD:
//                pos_x = parent_x + (10./14) * mult_x;
//                pos_y = parent_y + (10./14) * mult_y;
//                size = 18;
//                break;
            case SPIKE:
                pos_x = parent_x + 2 * mult_x;
                pos_y = parent_y + 2 * mult_y;
                size = 0;
                break;
            default:
                throw std::invalid_argument("hitbox can only be retrieved on CELL and SPIKE");
        }
    }
    pos_x = std::cos(this->rotation) * (pos_x - this->x) - std::sin(this->rotation) * (pos_y - this->y) + this->x;
    pos_y = std::sin(this->rotation) * (pos_x - this->x) + std::cos(this->rotation) * (pos_y - this->y) + this->y;
    return {pos_x, pos_y, size};
}

std::pair<double, double> User::get_direction_normalised(int component, User &to, int component2) const {
    auto [pos_x, pos_y, _] = this->get_hitbox(component);
    auto [pos2_x, pos2_y, size2] = to.get_hitbox(component2);
    double vec_x = pos2_x - pos_x, vec_y = pos2_y - pos_y;
    double vec_magnitude = sqrt(vec_x*vec_x + vec_y*vec_y);
    return {vec_x / vec_magnitude, vec_y / vec_magnitude};
}

bool User::collide_with_user(User& user) {
    if (distance(*this, user) > this->size + user.size) return false;

    return this->foreach_component<bool>([&](Component &c, int i, bool acc) -> bool {
        if (c.type != CELL && c.type != SPIKE) return acc;
        double pos_x, pos_y;
        int hitbox_size;
        std::tie(pos_x, pos_y, hitbox_size) = this->get_hitbox(i);
        return user.foreach_component<bool>([&](Component &c2, int j, bool acc2) -> bool {
            if (c2.type != CELL && c2.type != SPIKE) return acc2;
            auto [pos2_x, pos2_y, size2] = user.get_hitbox(j);
            if (distance(pos_x, pos_y, pos2_x, pos2_y) > hitbox_size + size2) {
                int layered_part = this->get_layered_part(i, pos2_x, pos2_y);
                int layered_part2 = user.get_layered_part(j, pos_x, pos_y);
                BODYPART_TYPE a = layered_part == -1 ? c.type : this->components[layered_part].type;
                BODYPART_TYPE b = layered_part2 == -1 ? c2.type : this->components[layered_part2].type;

                switch (a) {
                    case BOUNCE:
                        if (this->components[layered_part].is_working()) {
                            switch (b) {
                                case BOUNCE:
                                    if (this->components[layered_part].is_working()) {
                                        // TODO: tbd
                                        break;
                                    } // intentional fallthrough
                                case CELL:
                                case SHIELD:
                                    user.bounce_back(*this, layered_part);
                                    break;
                                case SPIKE:
                                    user.bounce_back(*this, layered_part, true);
                                    this->components[layered_part].inflated = 0;
                                    break;
                                case NONE:
                                    throw std::invalid_argument("bodypart with type NONE encountered during collision");
                            }
                            break;
                        } // intentional fallthrough
                    case CELL:
                        switch (b) {
                            case BOUNCE:
                                if (user.components[layered_part2].is_working()) {
                                    this->bounce_back(user, layered_part2);
                                    break;
                                } // intentional fallthrough
                            case CELL:
                            case SHIELD:
                                this->transfer_momentum(i, user, j);
                                user.transfer_momentum(j, *this, i);
                                break;
                            case SPIKE:
                                if (this->damage(i, SPIKE_DMG)) {
                                    user.kills++;
                                }
                                break;
                            case NONE:
                                throw std::invalid_argument("bodypart with type NONE encountered during collision");
                        }
                        break;
                    case SPIKE:
                        switch (b) {
                            case BOUNCE:
                                if (user.components[layered_part2].is_working()) {
                                    this->bounce_back(user, layered_part2, true);
                                    user.components[layered_part2].inflated = 0;
                                    break;
                                } // intentional fallthrough
                            case CELL:
                                if (user.damage(j, SPIKE_DMG)) {
                                    this->kills++;
                                }
                                break;
                            case SPIKE:
                                if (user.damage(j, SPIKE_DMG)) {
                                    this->kills++;
                                }
                                if (this->damage(i, SPIKE_DMG)) {
                                    user.kills++;
                                }
                                break;
                            case SHIELD:
                                this->transfer_momentum(i, user, j);
                                user.transfer_momentum(j, *this, i);
                                break;
                            case NONE:
                                throw std::invalid_argument("bodypart with type NONE encountered during collision");
                        }
                        break;
                    case SHIELD:
                        switch (b) {
                            case BOUNCE:
                                if (user.components[layered_part2].is_working()) {
                                    this->bounce_back(user, layered_part2);
                                    break;
                                } // intentional fallthrough
                            case CELL:
                            case SPIKE:
                            case SHIELD:
                                this->transfer_momentum(i, user, j);
                                user.transfer_momentum(j, *this, i);
                                break;
                            case NONE:
                                throw std::invalid_argument("bodypart with type NONE encountered during collision");
                        }
                        break;
                    case NONE:
                        throw std::invalid_argument("bodypart with type NONE encountered during collision");
                }
                return true;
            }
            return acc2;
        }, acc);

    }, false);
}

void User::bounce_back(User &bouncer, int bounce_component, bool popped) {
    // we have an off-by-one error here on purpose.
    // one of the faces _has_ to contain the bounce. If it gets to the 7th face and throws
    // an exception due to the illegal array access something has already gone wrong
    int face = 0;
    for (; face <= 7; face++) {
        if (bouncer.components[bouncer.components[bounce_component].body].faces[face] == bounce_component) break;
    }
    double x_magn = std::cos(FACE_TO_ANGLE[face]), y_magn = std::sin(FACE_TO_ANGLE[face]);
    this->momentum_x += BOUNCE_MOMENTUM * (popped ? BOUNCE_POPPED_MULTIPLIER : 1) * x_magn;
    this->momentum_y += BOUNCE_MOMENTUM * (popped ? BOUNCE_POPPED_MULTIPLIER : 1) * y_magn;
    bouncer.momentum_x += BOUNCE_MOMENTUM_SELF * (popped ? BOUNCE_POPPED_MULTIPLIER : 1) * -x_magn;
    bouncer.momentum_y += BOUNCE_MOMENTUM_SELF * (popped ? BOUNCE_POPPED_MULTIPLIER : 1) * -y_magn;
}

void User::transfer_momentum(int component, User &user, int other_component) {
    auto [pos_x, pos_y, _] = get_hitbox(component);
    auto [pos2_x, pos2_y, size2] = user.get_hitbox(other_component);
    double vec_x = pos2_x - pos_x, vec_y = pos2_y - pos_y;
    double vec_magnitude = sqrt(vec_x*vec_x + vec_y*vec_y);
    double vec_x_norm = vec_x / vec_magnitude, vec_y_norm = vec_y / vec_magnitude;

    if ((vec_x_norm < 0 && momentum_x < 0)
    || (vec_x_norm >= 0 && momentum_x >= 0)) {
        double amt_to_transfer = momentum_x * MOMENTUM_TRANSFER_ON_COLLISION * vec_x_norm;
        user.momentum_x += amt_to_transfer;
        this->momentum_x -= amt_to_transfer;
    }
    if ((vec_y_norm < 0 && momentum_y < 0)
    || (vec_y_norm >= 0 && momentum_y >= 0)) {
        double amt_to_transfer = momentum_y * MOMENTUM_TRANSFER_ON_COLLISION * vec_y_norm;
        user.momentum_y += amt_to_transfer;
        this->momentum_y -= amt_to_transfer;
    }
}

int User::get_layered_part(int component, double to_x, double to_y) const {
    if (this->components[component].type != CELL) return NONE;
    auto [pos_x, pos_y, _] = this->get_hitbox(component);

    double vec_x = to_x - pos_x, vec_y = to_y - pos_y;

    double angle = std::atan2(vec_y, vec_x);
    double local_angle = this->rotation + angle;

    // this (implicit) cast is safe - after flooring we're guaranteed to have an integer number which should
    // convert to int without issues
    int offset = std::floor(local_angle * 3 / M_PI) + 3; // NOLINT(cppcoreguidelines-narrowing-conversions)
    if (offset == 6) offset = 0;

    return this->components[component].faces[offset];
}

bool User::collide_with_resource(Resource &res) {
    if (res.amt <= 0) return false;
    double resource_size = std::min((double) RESOURCE_SIZE_MAX, RESOURCE_SIZE_MIN + res.amt * RESOURCE_SIZE_STEP);
    if (distance(this->x, this->y, res.x, res.y) > this->size + resource_size) return false;

    return this->foreach_component<bool>([&](const Component &c, int i, bool acc) -> bool {
        if (c.type != CELL && c.type != SPIKE) return acc;

        auto [pos_x, pos_y, hitbox_size] = this->get_hitbox(i);
        if (distance(pos_x, pos_y, res.x, res.y) > resource_size + hitbox_size) return acc;

        int layered_component = this->get_layered_part(i, res.x, res.y);
        if (layered_component != -1 && this->components[layered_component].type == SHIELD) return acc;

        int amt = std::min(res.amt, MINING_RATE * (c.type == SPIKE ? SPIKE_MINING_MULTIPLIER : 1));
        this->resources += amt;
        res.amt -= amt;
        return true;
    }, false);
}

/**
 * Returns an approximation of the user's size as measured by distance between root and component farthest from it.
 * @return an overestimate (by no more than 2*INNER_CELL_RADIUS) of the distance between the root cell and the outermost
 *         edge of the outermost bodypart
 */
double User::get_size() const {
    return size;
}

void User::compute_size() {
    this->size = this->foreach_component<double>([this](Component &c, double dist) -> double {
        if (c.type != CELL) return dist;
        auto [pos_x, pos_y] = rel_pos(this->x, this->y, std::get<0>(c.coords),
                                      std::get<1>(c.coords), std::get<2>(c.coords));
        double new_distance = distance(pos_x, pos_y, this->x, this->y);
        return std::max(new_distance, dist);
        }, 0.);
}

template<typename T> T User::foreach_component(const std::function<T(Component&, T)>& f, T init) {
    return this->foreach_component<T>([&](Component&c, int i, T a) {
        return f(c, a);
    }, init);
}
void User::foreach_component(const std::function<void(Component&)>& f) const{
    this->foreach_component([&](Component &c, int) {
        f(c);
    });
}
template<typename T> T User::foreach_component(const std::function<T(Component&, int, T)>& f, T init) {
    T acc = init;
    for (int i = 0; i < this->next_component; i++) {
        acc = f(this->components[i], i, acc);
    }
    return acc;
}
void User::foreach_component(const std::function<void(Component&, int)>& f) const{
    for (int i = 0; i < this->next_component; i++) {
        f(this->components[i], i);
    }
}

