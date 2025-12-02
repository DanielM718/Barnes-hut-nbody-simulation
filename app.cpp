#include <stdio.h>
#include "node.h"

int main(){
    node* root = new node(100);
    root->AddObject(1, vector(0,0,0), vector(0,0,0));
    root->AddObject(0.01, vector(1,0,0), vector(0, 2*M_PI, 0));

    while (1){
        root->simulate(root);
        root->rebuild(root);
        root->animate(root);
        printf("F\n");
    }
    
    return 1;
}