#ifndef HEXARENA_BACKEND_COMPONENT_H
#define HEXARENA_BACKEND_COMPONENT_H

#include "constants.h"

struct Component{
    BODYPART_TYPE type;
    /// order is: top-right, top, top-left, bottom-left, bottom, bottom-right. Counterclockwise starting at the x axis
    int faces[6];
    int health;
    bool dead;
    std::tuple<int, int, int> coords;

    // bounce
    int inflated;

    [[nodiscard]] bool is_working() const {
        return inflated >= MAX_INFLATE;
    }

    // bounce/spike/shield
    int body;
};

#endif //HEXARENA_BACKEND_COMPONENT_H
