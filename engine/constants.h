#ifndef HEXARENA_BACKEND_CONSTANTS_H
#define HEXARENA_BACKEND_CONSTANTS_H

#include <map>

enum BODYPART_TYPE {
    CELL,
    SPIKE,
    SHIELD,
    BOUNCE
};

std::map<BODYPART_TYPE, int> BODYPART_COST = {{CELL, 20},
                                                     {SPIKE, 15},
                                                     {SHIELD, 10},
                                                     {BOUNCE, 15}};

#define MAX_HEALTH 100
#define MAX_INFLATE 100
#define INFLATE_RATE 2
#define REGEN_RATE 1
#define SPIKE_DMG 3
#define BANANA true

enum ACTION {
    MOVE,
    DESTROY,
    COLLIDE
};

#define CHEATS_ENABLED 1

#define RESOURCE_SIZE_STEP 10
#define RESOURCE_SIZE_MIN 5
#define RESOURCE_SIZE_MAX 20

#define MINING_RATE 1

#endif //HEXARENA_BACKEND_CONSTANTS_H
