#include <stdio.h>
#include "node.h"

int main(){
    node* root = new node(100);
    root->AddObject(0.01, vector(0,0,0), vector(0,0,0));
    root->AddObject(0.01, vector(0.5,0,0), vector(0,0,0));

    while (1){
        root->simulate(root);
        root->rebuild(root);
        root->animate(root);
        //printf("center3 0 0 0\n");
        printf("F\n");
    }
    
    return 1;
}