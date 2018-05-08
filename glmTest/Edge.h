//
//  Edge.hpp
//  glmTest
//
//  Created by Lucas Andrade on 30/01/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#ifndef Edge_h
#define Edge_h

#include <stdio.h>

#include "Vector.h"

namespace primitives {
    class Edge{
    public:
        Edge();
        Edge(Vector *v0, Vector *v1);
        ~Edge();
        
        Vector *v0, *v1;

        Face *leftFace;
        Face *rightFace;
        
        bool toBeErased;
        Vector *getVector(int i);
    };
}

#endif /* Edge_h */
