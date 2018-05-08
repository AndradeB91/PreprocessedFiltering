//
//  Boundary.cpp
//  glmTest
//
//  Created by Lucas Andrade on 01/02/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#include "Boundary.h"

using namespace Data;

Boundary::Boundary(){
    
}

Boundary::~Boundary(){
    
}

std::list<Vector> Boundary::getOutlineVectors(){
    for(std::list<Edge>::iterator it = outline.begin(); it != outline.end(); it++){
        outlineVectors.push_back( *(*it).v0 );
    }
    
    return this->outlineVectors;
}
