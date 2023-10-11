#pragma once

#define EPS 1e-9

template <class T>
bool rough_eq(T lhs, T rhs, T epsilon = EPS) { 
    return fabs(lhs - rhs) < epsilon;
}

template <class T>
bool rough_lt(T lhs, T rhs, T epsilon = EPS) {
    return rhs - lhs >= epsilon;
}

template <class T>
bool rough_lte(T lhs, T rhs, T epsilon = EPS) {
    return rhs - lhs > -epsilon;
}