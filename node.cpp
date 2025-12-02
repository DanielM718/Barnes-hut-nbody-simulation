#include "node.h"
#include <math.h>
#include <stdio.h>

node::node(double width) : 
    space(width),
    one(nullptr),
    two(nullptr),
    three(nullptr),
    four(nullptr),
    five(nullptr),
    six(nullptr),
    seven(nullptr),
    eight(nullptr),
    state(leaf),
    id(-1)
    {
}

node::node(double width, int objects, double mass, vector com, vector center):
    space(width, objects, mass, com, center),
    one(nullptr),
    two(nullptr),
    three(nullptr),
    four(nullptr),
    five(nullptr),
    six(nullptr),
    seven(nullptr),
    eight(nullptr),
    state(leaf),
    id(-1)
    {
}

node::node(double width, int objects, double mass, vector com, vector center, int id):
    space(width, objects, mass, com, center),
    one(nullptr),
    two(nullptr),
    three(nullptr),
    four(nullptr),
    five(nullptr),
    six(nullptr),
    seven(nullptr),
    eight(nullptr),
    state(leaf),
    id(id)
    {
}

node::~node(){
    delete one;
    delete two;
    delete three;
    delete four;
    delete five;
    delete six;
    delete seven;
    delete eight;
}

void node::SetOne(node* one){this->one = one;}
void node::SetTwo(node* two){this->two = two;}
void node::SetThree(node* three){this->three = three;}
void node::SetFour(node* four){this->four = four;}
void node::SetFive(node* five){this->five = five;}
void node::SetSix(node* six){this->six = six;};
void node::SetSeven(node* seven){this->seven = seven;}
void node::SetEight(node* eight){this->eight = eight;}

const double node::G = 4*(M_PI*M_PI);
const double node::dt = 0.0005;
const double node::theta = 0.5;
Body* node::leafNodes[100000] = {nullptr};
int node::memberCount = 0;


// for the Barnes Hut method I am going to construct a tree
// note dynamically changing the tree as the objects move in space is a lot of work
// Maybe 
void node::AddObject(double mass, vector com, vector v){
    if((state == leaf) && (objects == 0)){
        // if its empty this will be the first body
        this->id = memberCount;
        leafNodes[id] = new Body(mass, com, v);
        memberCount++;

        SetObject(1);
        SetMass(mass);
        SetCom(com);

        return; // nothing else needs to be done
    }

    if((state == leaf) && (objects == 1)){
        // if it has more than one object its an internal node now
        state = internal;

        // reset this internal node to a blank slate
        SetMass(0.0);
        SetCom(vector(0, 0, 0));
        SetObject(0);

        InsertChild(leafNodes[id]->mass, leafNodes[id]->position, leafNodes[id]->velocity, id);
        id = -1;
    }
    InsertChild(mass, com, v, id);

    // I want to avoid updating these values recursively because theyre very annoying to debug and keep track of
    // this way is more straight forward.
    double totalMass = 0.0;
    vector weightedSum(0, 0, 0);
    int totalObjects = 0;

    // after inserting the new node I can add up all the masses again.
    node* kids[8] = { one, two, three, four, five, six, seven, eight };

    for (int i = 0; i < 8; i++){
        node* c = kids[i];
        if((c != nullptr) && (c->mass > 0.0)){
            totalMass += c->mass;
            weightedSum += c->com * c->mass;
            totalObjects += c->objects;
        }
    }
    SetCom(weightedSum / totalMass);
    SetMass(totalMass);
    SetObject(totalObjects);
}

void node::AddExisting(int id){
    Body* b = leafNodes[id];
    vector com = b->position;
    vector v = b->velocity;
    double mass = b->mass;

    if((state == leaf) && (objects == 0)){
        SetCom(com);
        SetMass(mass);
        SetObject(1);
        this->id = id;
        return;
    }

    if((state == leaf) && (objects == 1)){
        // if it has more than one object its an internal node now
        state = internal;

        double m_old = this->mass;
        vector com_old = this->com;
        vector velocity = leafNodes[this->id]->velocity;

        // reset this internal node to a blank slate
        SetMass(0.0);
        SetCom(vector(0, 0, 0));
        SetObject(0);


        InsertExisting(m_old, com_old, velocity, this->id);
        this->id = -1;
    }
    InsertExisting(mass, com, v, id);

    // I want to avoid updating these values recursively because theyre very annoying to debug and keep track of
    // this way is more straight forward.
    double totalMass = 0.0;
    vector weightedSum(0, 0, 0);
    int totalObjects = 0;

    // after inserting the new node I can add up all the masses again.
    node* kids[8] = { one, two, three, four, five, six, seven, eight };

    for (int i = 0; i < 8; i++){
        node* c = kids[i];
        if((c != nullptr) && (c->mass > 0.0)){
            totalMass += c->mass;
            weightedSum += c->com * c->mass;
            totalObjects += c->objects;
        }
    }
    SetCom(weightedSum / totalMass);
    SetMass(totalMass);
    SetObject(totalObjects);
    
}

int node::SetID(double m, vector r, vector v){
    if(id == -1){
        int tmp = memberCount;
        leafNodes[memberCount] = new Body(m, r, v);
        memberCount++;
        return tmp;
    }
    return id;
}

void node::InsertChild(double m, vector r, vector v, int id){
    double x = this->center.x;
    double y = this->center.y;
    double z = this->center.z;

    double cx, cy, cz;


    // child width and subdivison offset
    double childWidth = width / 2.0;
    double offset = childWidth / 2.0;


    
    // this handles which branch the new object will be in the oct tree
    if(r.x >= x && r.y >= y && r.z >= z){
        // top, up, right of cube
        // if the branch is empty I can just append it to the end
        if(one == nullptr){
            cx = x + offset;
            cy = y + offset;
            cz = z + offset;
            if(id == -1){
                id = SetID(m, r, v);
            }
            SetOne(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        } // but if there is another node, then that means it must traverse to another node recursively
        else{one->AddObject(m, r, v);}
    }
    else if (r.x < x && r.y >= y && r.z >= z){
        // top, up, left of cube
        if(two == nullptr){
            cx = x - offset;
            cy = y + offset;
            cz = z + offset;
            if(id == -1){
                id = SetID(m, r, v);
            }
            SetTwo(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{two->AddObject(m, r, v);}
    }
    else if(r.x >= x && r.y < y && r.z >= z){
        // top, down, right of cube
        if(three == nullptr){
            cx = x + offset;
            cy = y - offset;
            cz = z + offset;
            if(id == -1){
                id = SetID(m, r, v);
            }
            SetThree(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{three->AddObject(m, r, v);}
    }
    else if(r.x < x && r.y < y && r.z >= z){
        // top, down, left of cube
        if(four == nullptr){
            cx = x - offset;
            cy = y - offset;
            cz = z + offset;
            if(id == -1){
                id = SetID(m, r, v);
            }
            SetFour(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{four->AddObject(m, r, v);}
    }
    else if(r.x >= x && r.y >= y && r.z < z){
        // bottom, up, right of cube
        if(five == nullptr){
            cx = x + offset;
            cy = y + offset;
            cz = z - offset;
            if(id == -1){
                id = SetID(m, r, v);
            }
            SetFive(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{five->AddObject(m, r, v);}
    }
    else if(r.x < x && r.y >= y && r.z < z){
        // bottom, up, left of cube
        if(six == nullptr){
            cx = x - offset;
            cy = y + offset;
            cz = z - offset;
            if(id == -1){
                id = SetID(m, r, v);
            }
            SetSix(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{six->AddObject(m, r, v);}
    }
    else if(r.x >= x && r.y < y && r.z < z){
        // bottom, down, right of cube
        if(seven == nullptr){
            cx = x + offset;
            cy = y - offset;
            cz = z - offset;
            if(id == -1){
                id = SetID(m, r, v);
            }
            SetSeven(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{seven->AddObject(m, r, v);}
    }
    else{
        // bottom, down, left of cube
        if(eight == nullptr){
            cx = x - offset;
            cy = y - offset;
            cz = z - offset;
            if(id == -1){
                id = SetID(m, r, v);
            }
            SetEight(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{eight->AddObject(m, r, v);}
    }
}

void node::InsertExisting(double m, vector r, vector v, int id){
    double x = this->center.x;
    double y = this->center.y;
    double z = this->center.z;

    double cx, cy, cz;


    // child width and subdivison offset
    double childWidth = width / 2.0;
    double offset = childWidth / 2.0;


    
    // this handles which branch the new object will be in the oct tree
    if(r.x >= x && r.y >= y && r.z >= z){
        // top, up, right of cube
        // if the branch is empty I can just append it to the end
        if(one == nullptr){
            cx = x + offset;
            cy = y + offset;
            cz = z + offset;
            SetOne(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        } // but if there is another node, then that means it must traverse to another node recursively
        else{one->AddExisting(id);}
    }
    else if (r.x < x && r.y >= y && r.z >= z){
        // top, up, left of cube
        if(two == nullptr){
            cx = x - offset;
            cy = y + offset;
            cz = z + offset;
            SetTwo(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{two->AddExisting(id);}
    }
    else if(r.x >= x && r.y < y && r.z >= z){
        // top, down, right of cube
        if(three == nullptr){
            cx = x + offset;
            cy = y - offset;
            cz = z + offset;
            SetThree(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{three->AddExisting(id);}
    }
    else if(r.x < x && r.y < y && r.z >= z){
        // top, down, left of cube
        if(four == nullptr){
            cx = x - offset;
            cy = y - offset;
            cz = z + offset;
            SetFour(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{four->AddExisting(id);}
    }
    else if(r.x >= x && r.y >= y && r.z < z){
        // bottom, up, right of cube
        if(five == nullptr){
            cx = x + offset;
            cy = y + offset;
            cz = z - offset;
            SetFive(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{five->AddExisting(id);}
    }
    else if(r.x < x && r.y >= y && r.z < z){
        // bottom, up, left of cube
        if(six == nullptr){
            cx = x - offset;
            cy = y + offset;
            cz = z - offset;
            SetSix(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{six->AddExisting(id);}
    }
    else if(r.x >= x && r.y < y && r.z < z){
        // bottom, down, right of cube
        if(seven == nullptr){
            cx = x + offset;
            cy = y - offset;
            cz = z - offset;
            SetSeven(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{seven->AddExisting(id);}
    }
    else{
        // bottom, down, left of cube
        if(eight == nullptr){
            cx = x - offset;
            cy = y - offset;
            cz = z - offset;
            SetEight(new node(childWidth, 1, m, r, vector(cx, cy, cz), id));
        }
        else{eight->AddExisting(id);}
    }    
}

// misleading as I compute the acceleration only but it saves me some time
// there is one optimization i will implement later? 
void node::computeForce(const node* n, const vector& com, const double theta, vector& alpha){
    if(!n) {return;} // since Ill be doing this recursively, it will stop
    // when it reaches a leaf node 

    double s = n->width; // in the Barnes Hut method the approximation is controlled by s/d
    // the

    vector r = n->com - com; // compute the displacment between the two nodes center of mass
    double d = this->magnitude(r);

    if(n->objects == 1){ // if its one then it means we are at a leaf node
        // in that case we can just direct force and return
        if(d == 0) {return;} // if d is 0 then it means we are at the body we are comparing so skip
        alpha += ((G*(n->mass)) / (d*d*d))*r;
        return; 
    }

    // now this is the logic of the actual implementation of the Barnes Hut method
    // if the angle is small enough then we can treat it as one object
    // how high theta is determines the accuracy of the simualtion
    // if its 0, then it will compute it as sum all
    if(s/d < theta){
        if(d == 0) {return;}
        alpha += ((G*(n->mass)) / (d*d*d))*r;
        return;
    }
    else{
        // now recursion time.
        computeForce(n->one, com, theta, alpha);
        computeForce(n->two, com, theta, alpha);
        computeForce(n->three, com, theta, alpha);
        computeForce(n->four, com, theta, alpha);
        computeForce(n->five, com, theta, alpha);
        computeForce(n->six, com, theta, alpha);
        computeForce(n->seven, com, theta, alpha);
        computeForce(n->eight, com, theta, alpha);
    }
}


double node::magnitude(vector tmp){
    return sqrt(tmp.x*tmp.x + tmp.y*tmp.y + tmp.z*tmp.z);
}

void node::PositionHalfStep(){
    leafNodes[id]->position += leafNodes[id]->velocity*dt/2;
}

void node::VelocityHalfStep(){
    
    leafNodes[id]->velocity += leafNodes[id]->alpha*dt;
}

void node::traversal(node* root, int opp){
    // start at root so first iteration root=curr
    // must be a leaf and implicitly also have 1 object
    if((state == leaf) && (objects == 1)){
        // when it is a leaf then we simply compute the force
        switch (opp)
        {
        case forceV:
            leafNodes[id]->alpha = vector(0, 0, 0);
            computeForce(root, leafNodes[id]->position, theta, leafNodes[id]->alpha);
            VelocityHalfStep();
            break;
        case Pstep:
            PositionHalfStep();
            break;
        default:
            printf("bad opp");
            break;
        }
        return;
    }
    if(this->one){// its important I keep track of a pointer to the root
        // it will essentially travel throughout the tree until it hits every ndoe
        // if its an internal node it will visit each child node ones
        // it does this by rerunning simulate and passing it the pointer to the next node
        this->one->traversal(root, opp);
    }
    if(this->two){
        this->two->traversal(root, opp);
    }
    if(this->three){
        this->three->traversal(root, opp);
    }
    if(this->four){
        this->four->traversal(root, opp);
    }
    if(this->five){
        this->five->traversal(root, opp);
    }
    if(this->six){
        this->six->traversal(root, opp);
    }
    if(this->seven){
        this->seven->traversal(root, opp);
    }
    if(this->eight){
        this->eight->traversal(root, opp);
    }
}

void node::simulate(node* root){
    // half step
    root->traversal(root, Pstep);

    // velocity step
    root->traversal(root, forceV);

    // second half step
    root->traversal(root, Pstep);

}

void node::rebuild(node*& root){
    if(!root){return;} // if root pointer is a null pointer stop

    double width = root->width;
    vector center = root->center;

    node* newRoot = new node(width);
    newRoot->SetCenter(center);

    for (int i = 0; i < memberCount; i++){
        if(leafNodes[i]){
            newRoot->AddExisting(i);
        }
    }

    delete root;
    root = newRoot;
}


