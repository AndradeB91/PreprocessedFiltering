//
//  Vector.cpp
//  glmTest
//
//  Created by Lucas Andrade on 17/10/16.
//  Copyright © 2016 Lucas Andrade. All rights reserved.
//

#include "Vector.h"

using namespace primitives;

Vector::Vector(){
    toBeErased = false;
}

Vector::Vector(float x, float y, float z){
    toBeErased = false;
    this->x = x;
    this->y = y;
    this->z = z;
}

Vector::~Vector(){}

double Vector::distance(Vector v){
    double exp = pow(this->x - v.x, 2) + pow(this->y - v.y, 2) + pow(this->z - v.z, 2);
    return sqrt(exp);
}

float Vector::getCoord(int i){
    if(i == 0) return this->x;
    else if(i == 1) return this->y;
    else if(i == 2) return this->z;
    else return -1;
}

Vector Vector::normalize(){
    Vector v;
    float den = sqrt((this->x * this->x) + (this->y * this->y) + (this->z * this->z) );
    v.x = this->x/den;
    v.y = this->y/den;
    v.z = this->z/den;
    return v;
}

float Vector::length(){
    return sqrt((this->x * this->x) + (this->y * this->y) + (this->z * this->z) );
}

float Vector::dot(Vector v){
    return ((this->x * v.x) + (this->y * v.y) + (this->z * v.z));
}
