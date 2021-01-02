#include "comp.h"


using namespace civ; 

int main() {
    
    printf("Program start:\n");
    
    frame f; 

    // element 1 
    f.add_element(new element(
        2e8,        // E
        0.00048,    // I 
        0.075,      // A
        5,          // L
        7880,       // rho
        0,          // x1
        0,          // y1
        0,          // x2
        5,          // y2
        90));       // theta

    // element 2 
    f.add_element(new element(
        2e8,        // E
        0.00048,    // I 
        0.075,      // A
        5,          // L
        7880,       // rho
        0,          // x1
        5,          // y1
        0,          // x2
        10,         // y2
        90));       // theta

    // element 3 
    f.add_element(new element(
        2e8,        // E
        0.00048,    // I 
        0.075,      // A
        6,          // L
        7880,       // rho
        0,          // x1
        5,          // y1
        6,          // x2
        5,          // y2
        0));        // theta

    f.global_k();

    printf("Program completed.\n");
}
