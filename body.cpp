#include "body.h"

Body::Body(double mass, vector com, vector velocity):
    mass(mass),
    position(com),
    velocity(velocity),
    alpha(vector(0, 0, 0)){}