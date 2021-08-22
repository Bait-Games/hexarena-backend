#ifndef HEXARENA_BACKEND_USER_H
#define HEXARENA_BACKEND_USER_H

#include <vector>
#include <set>

#include "constants.h"
#include "component.h"
#include "utils.h"

#define CELL_INNER_RADIUS 28

class User {
private:
    double size;
    /// Pointer to a function whose purpose is to drop the given amount of resources on the ground at the given position.
    std::function<void(int, double, double)> drop;

    void find_connected(int node, std::set<int> &set);

    template<typename T> T foreach_component(const std::function<T(Component&, T)>&, T init);
    void foreach_component(const std::function<void(Component&)>&) const;
    template<typename T> T foreach_component(const std::function<T(Component&, int , T)>&, T init);
    void foreach_component(const std::function<void(Component&, int)>&) const;

    [[nodiscard]] std::tuple<double, double, int> get_hitbox(int component) const;
    std::pair<double, double> get_direction_normalised(int component, User &to, int component2) const;
    [[nodiscard]] int get_layered_part(int component, double to_x, double to_y) const;
    void transfer_momentum(int component, User &user, int other_component);
    void bounce_back(User &bouncer, int bounce_component, bool popped = false);

    void compute_size();


public:
    int id;
    double x;
    double y;
    double momentum_x;
    double momentum_y;
    double rotation;
    Component* components;
    int components_size;
    int next_component;
    int resources;
    int kills;
    bool destroyed;
    bool moved;
    bool rotated;
    std::vector<void (*)(const User&)> observers;

#if CHEATS_ENABLED
    // char cheat_seq[20]; // TODO: cheats
#endif


    User(int id, double x, double y, std::function<void(int, double, double)> drop);
    void tick_reset();
    void tick_parts() const;
    int grow(int part, int face, BODYPART_TYPE type);
    bool damage(int part, int amt);
    bool shrink(int part);
    bool collide_with_user(User& user);
    bool collide_with_resource(Resource &res);
    double get_size() const;

};

#endif //HEXARENA_BACKEND_USER_H
