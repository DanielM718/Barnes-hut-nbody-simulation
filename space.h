#ifndef SPACE_H
#define SPACE_H
#include "vector.h"

class space{
    protected:
        double width; // width of cube
        int objects; // number of objects in space
        double mass; // mass of all objects
        vector com; // center of mass
        vector center; // center of space
    public:
        space(double width);
        space(double width, int objects, double mass, vector com, vector center);

        void SetWidth(double width);
        void SetObject(int objects);
        void IncObject();
        void SetMass(double mass);
        void AddMass(double mass);
        void SetCom(vector com);
        void SetCenter(vector center);
};

#endif