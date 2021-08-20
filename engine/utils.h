#ifndef HEXARENA_BACKEND_UTILS_H
#define HEXARENA_BACKEND_UTILS_H

#include "user.h"

struct Resource {
    double x;
    double y;
    int amt;
};

double distance(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
}

#endif //HEXARENA_BACKEND_UTILS_H
