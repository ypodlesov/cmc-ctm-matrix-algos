#pragma once

#include <cmath>

#define EPS 1e-9

template <class T>
bool RoughEq(T lhs, T rhs, T epsilon = EPS) { 
    return std::fabs(lhs - rhs) < epsilon;
}

template <class T>
bool RoughLT(T lhs, T rhs, T epsilon = EPS) {
    return rhs - lhs >= epsilon;
}

template <class T>
bool RoughtLTE(T lhs, T rhs, T epsilon = EPS) {
    return rhs - lhs > -epsilon;
}