#ifndef HEXARENA_BACKEND_ENGINE_H
#define HEXARENA_BACKEND_ENGINE_H

#include <random>
#include <thread>
#include "constants.h"
#include "user.h"
#include "utils.h"


class Engine {
private:
    std::random_device r;
    std::mt19937 gen{r()};

    long tick_num;
    std::map<int, User> users;
    int next_user_id;
    std::vector<Resource> resources;
    void tick();
    void run();
    int add_user(double x, double y);
    bool running;
    std::thread runner;
    std::mutex mutex;
    std::vector<void (*)(const std::vector<Resource>&, const std::map<int, User>&)> observers;

public:
    Engine();
    ~Engine();
    void move(int id, double angle);
    void rotate(int id, double angle);
    int create();
    void remove(int id);
    int attach(int id, BODYPART_TYPE type, int part, int face);
    int detach(int id, int part);
    const User & user(int id);
    [[deprecated("renamed to user()")]] const User & info(int id);
    void observe(void (*cb)(const std::vector<Resource>&, const std::map<int, User>&));
    [[deprecated("renamed to observe()")]]
    void register_global(void (*cb)(const std::vector<Resource>&, const std::map<int, User>&));
    void observe(int user, void (*cb)(const User &));
    void drop(int amt, double x, double y);
};

#endif //HEXARENA_BACKEND_ENGINE_H
