//
//  Face.cpp
//  glmTest
//
//  Created by Lucas Andrade on 30/09/16.
//  Copyright © 2016 Lucas Andrade. All rights reserved.
//

#include "Face.h"

using namespace primitives;

Face::Face(){
    this->toBeErased = false;
    this->processed = false;
}

Face::~Face(){}

float Face::area(){
    
    Vector v01 = *v1 - *v0;
    Vector v02 = *v2 - *v0;

    glm::vec3 glmv01, glmv02;
    glmv01.x = v01.x;
    glmv01.y = v01.y;
    glmv01.z = v01.z;
    
    glmv02.x = v02.x;
    glmv02.y = v02.y;
    glmv02.z = v02.z;

    glm::vec3 res = glm::cross(glmv01, glmv02);
    
    float r = glm::length(res);
    
    return r/2;
}

Vector Face::center(){
    Vector sum = *v0 + *v1 + *v2;
    sum.x = sum.x/3.0f;
    sum.y = sum.y/3.0f;
    sum.z = sum.z/3.0f;
    return sum;
}

Vector *Face::getVector(int i){
    if(i == 0) return this->v0;
    else if(i == 1) return this->v1;
    else return this->v2;
}

std::list<Face *> Face::getFacePatch(){
    
    std::list<Face *> patch;
    
    for(int i=0; i<3; i++){
        
        Vector *v = this->getVector(i);
        
        std::list<Face *>::iterator it;
        
        for(it = v->faceList.begin(); it != v->faceList.end(); it++){
            bool exists = false;
            
            if(patch.size() == 0){
                patch.push_back(*it);
            }
            else{
                std::list<Face *>::iterator itPatch;
                for(itPatch = patch.begin(); itPatch != patch.end(); itPatch++){
                    if( &(**it) == &(**itPatch) ){
                        exists = true;
                        break;
                    }
                    
                }
                if(exists == false)
                    patch.push_back(*it);
            }
        }
    }
    return patch;
}

std::list<Face *> Face::realFaceNeighbors(){
    
    std::list<Face *> patch = this->getFacePatch();
    std::list<Face *> realNeighbors;
    std::list<Face *>::iterator i;
    
    for(i = patch.begin(); i != patch.end(); i++){
        int cont = 0;
        
        if( (**i).v0 == this->v0 || (**i).v1 == this->v0 || (**i).v2 == this->v0){
            cont++;
        }
        
        if( (**i).v0 == this->v1 || (**i).v1 == this->v1 || (**i).v2 == this->v1){
            cont++;
        }
        
        if( (**i).v0 == this->v2 || (**i).v1 == this->v2 || (**i).v2 == this->v2){
            cont++;
        }
        
        //It found a real Face neighbor
        if(cont == 2) {
            realNeighbors.push_back(&(**i));
        }
    }
    return realNeighbors;
}

float Face::signedVolume(){
    Vector cross = *this->v1 % *this->v2;
    return ((*this->v0).dot(cross))/6.0f;
}

float Face::quality(){
    float a = (*this->v1 - *this->v0).length();
    float b = (*this->v2 - *this->v1).length();
    float c = (*this->v0 - *this->v2).length();
    float k = 0.5*(a+b+c);
    
    float outRadius = (a*b*c)/(float)sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c));
    float inRadius  = (float)sqrt(k*(k-a)*(k-b)*(k-c))/k;
    
    return outRadius/inRadius;
}


void Face::print(){
    std::cout << this->v0->x << " " << this->v0->y << " " << this->v0->z << std::endl;
    std::cout << this->v1->x << " " << this->v1->y << " " << this->v1->z << std::endl;
    std::cout << this->v2->x << " " << this->v2->y << " " << this->v2->z << "\n" << std::endl;
}












