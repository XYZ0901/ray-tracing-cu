//
// Created by root on 1/29/23.
//

#ifndef CU_RTWEEKEND_H
#define CU_RTWEEKEND_H

#include <cmath>
#include <limits>
#include <memory>

using std::shared_ptr;
using std::make_shared;
using std::sqrt;

const float infinity = std::numeric_limits<float>::infinity();
const float pi = 3.14159265358979;

__device__ inline float degrees_to_radians(float degrees){
    return degrees * pi / 180.0f;
}

inline float clamp(float x,float min,float max){
    if(x < min) return min;
    if (x > max) return max;
    return x;
}

#include "vec3.h"
#include "ray.h"

#endif //CU_RTWEEKEND_H
