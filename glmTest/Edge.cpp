//
//  Edge.cpp
//  glmTest
//
//  Created by Lucas Andrade on 30/01/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#include "Edge.h"

using namespace primitives;

Edge::Edge(){
    this->leftFace = NULL;
    this->rightFace = NULL;
    this->toBeErased = false;
}

Edge::~Edge(){
    
}

Edge::Edge(Vector *v0, Vector *v1){
    this->v0 = v0;
    this->v1 = v1;
    this->leftFace = NULL;
    this->rightFace = NULL;
}

Vector *Edge::getVector(int i){
    if(i == 0) return this->v0;
    else if(i == 1) return this->v1;
    else return NULL;
}
