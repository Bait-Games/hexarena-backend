#ifndef HEXARENA_BACKEND_CONSTANTS_H
#define HEXARENA_BACKEND_CONSTANTS_H

#include <map>
#include <cmath>

enum BODYPART_TYPE {
    NONE = 0,
    CELL,
    SPIKE,
    SHIELD,
    BOUNCE
};

const std::map<BODYPART_TYPE, int> BODYPART_COST = {{CELL, 20},
                                                     {SPIKE, 15},
                                                     {SHIELD, 10},
                                                     {BOUNCE, 15}};
const double FACE_TO_ANGLE[6] = {M_PI/6, M_PI*3/6, M_PI*5/6, M_PI*7/6, M_PI*9/6, M_PI*11/6};

#define MAX_INFLATE 100
#define INFLATE_RATE 2
#define BOUNCE_MOMENTUM 3
#define BOUNCE_MOMENTUM_SELF 0
#define BOUNCE_POPPED_MULTIPLIER 2
#define MAX_HEALTH 100
#define REGEN_RATE 1
#define SPIKE_DMG 3

#define CHEATS_ENABLED 1

#define RESOURCE_SIZE_STEP 0.1
#define RESOURCE_SIZE_MIN 5
#define RESOURCE_SIZE_MAX 20

#define MINING_RATE 1
#define SPIKE_MINING_MULTIPLIER 2

#define MOMENTUM_GAIN_PER_TICK 1
#define MOMENTUM_MULTIPLIER_PER_TICK 0.8
/// Maximum amount we can move an entity without checking for collisions. Collisions are very expensive so we want
/// to increase this number as much as possible.
#define MAX_SAFE_MOVEMENT 2

#define MOMENTUM_TRANSFER_ON_COLLISION 0.5

#define TICK_RATE 30
#define WORLD_WIDTH 3000
#define WORLD_HEIGHT 3000
#define RESOURCE_DENSITY 3
#define NATURAL_RESOURCE_AMOUNT 5
#define MAX_ROTATION_ITERATIONS 100

#endif //HEXARENA_BACKEND_CONSTANTS_H
