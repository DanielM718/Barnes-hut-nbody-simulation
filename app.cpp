#include <stdio.h>
#include "node.h"

int main(){
    node rootNode(100); // defines simulation region 100 light years wide center at 0
    rootNode.AddObject(1, vector(0,0,0));
    rootNode.AddObject(1, vector(1, 0, 0));
    
}