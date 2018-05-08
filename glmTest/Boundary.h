//
//  Boundary.h
//  glmTest
//
//  Created by Lucas Andrade on 01/02/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#ifndef Boundary_h
#define Boundary_h

#include <stdio.h>

#include "Edge.h"

using namespace primitives;

namespace Data {
    class Boundary{
    public:
        Boundary();
        ~Boundary();
        
        std::list<Edge> outline;
        
        std::list<Vector> getOutlineVectors();
                
    private:
        
        
        std::list<Vector> outlineVectors;
    };
}

#endif /* Boundary_h */
