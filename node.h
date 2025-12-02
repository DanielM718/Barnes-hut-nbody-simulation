#ifndef NODE_H
#define NODE_H
#include "space.h"
#include "body.h"

// this will be how I differentiate leafs vs internal nodes
enum state_{
    internal = 0,
    leaf = 1
};

enum opp_{
    forceV = 0,
    Pstep = 1,
    anim = 2
};



class node: public space{
    private:
        static Body* leafNodes[100000];
        static int memberCount;
        node* one;
        node* two;
        node* three;
        node* four;
        node* five;
        node* six;
        node* seven;
        node* eight;

        static const double G;
        static const double dt;
        static const double theta;


        // COM will serve as a position in space class


        int state;

        int id;
    public:
        node(double width);
        node(double width, int objects, double mass, vector com, vector center);
        node(double width, int objects, double mass, vector com, vector center, int id);
        ~node();

        void SetOne(node* one);
        void SetTwo(node* two);
        void SetThree(node* three);
        void SetFour(node* four);
        void SetFive(node* five);
        void SetSix(node* six);
        void SetSeven(node* seven);
        void SetEight(node* eight);

        void AddObject(double mass, vector com, vector v);
        void AddExisting(int id);

        void InsertChild(double m, vector r, vector v, int id);
        void InsertExisting(double m, vector r, vector v, int id);

        void computeForce(const node* n, const vector& com, double theta, vector& alpha);
        void PositionHalfStep();
        void VelocityHalfStep();

        double magnitude(vector tmp);

        void simulate(node* root);
        void traversal(node* root, int opp);
        void rebuild(node*& root);
        void animate(node* root);

        void printCoords();
        void trail();

        int SetID(double m, vector r, vector v);
};

#endif