#include "space.h"

space::space(double width): width(width), objects(0), mass(0), com(vector(0, 0, 0)), center(vector(0, 0, 0)){}
space::space(double width, int objects, double mass, vector com, vector center): width(width), objects(objects), mass(mass), com(com), center(center){}

void space::SetWidth(double width){this->width = width;}
void space::SetObject(int objects){this->objects = objects;}
void space::IncObject(){this->objects++;};
void space::SetMass(double mass){this->mass = mass;}
void space::AddMass(double mass){this->mass += mass;}
void space::SetCom(vector com){this->com = com;}
void space::SetCenter(vector center){this->center = center;}