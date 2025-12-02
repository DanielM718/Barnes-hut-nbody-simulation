#ifndef BODY_H
#define BODY_H
#include "vector.h"

// so i figured this might be easier
// this allows me to seperate the 
struct Body {
    vector alpha;
    vector velocity;
    vector position;
    double mass;
    Body(double mass, vector com, vector velocity);
};
#endif