//
//  Vector.hpp
//  glmTest
//
//  Created by Lucas Andrade on 17/10/16.
//  Copyright © 2016 Lucas Andrade. All rights reserved.
//

#ifndef Vector_h
#define Vector_h

#include <stdio.h>
#include <list>

#include <math.h>

namespace primitives {
    
    class Face;
    
    class Vector{
    public:
        Vector();
        Vector(float x, float y, float z);
        ~Vector();
        
        float x,y,z;
        int id;
        
        std::list<Face *> faceList;
        
        //operators
        Vector operator-(const Vector v) {
            Vector vector;
            vector.x = this->x - v.x;
            vector.y = this->y - v.y;
            vector.z = this->z - v.z;
            return vector;
        }
        
        Vector operator+(const Vector v) {
            Vector vector;
            vector.x = this->x + v.x;
            vector.y = this->y + v.y;
            vector.z = this->z + v.z;
            return vector;
        }
//        
//        Vector operator*(float f){
//            Vector vector;
//            vector.x *= f;
//            vector.y *= f;
//            vector.z *= f;
//            return vector;
//        }
//        
//        Vector operator/(float f){
//            Vector vector;
//            vector.x /= f;
//            vector.y /= f;
//            vector.z /= f;
//            return vector;
//        }
//        
//        Vector operator/(int f){
//            Vector vector;
//            vector.x /= f;
//            vector.y /= f;
//            vector.z /= f;
//            return vector;
//        }
        
        //cross product
        Vector operator%(const Vector v){
            Vector vector;
            float uvi = this->y * v.z - v.y * this->z;
            float uvj = v.x * this->z - this->x * v.z;
            float uvk = this->x * v.y - v.x * this->y;
            vector.x = uvi;
            vector.y = uvj;
            vector.z = uvk;
            return vector;
        }
    
        
        Vector normalize();
        
        float length();
        
        float dot(Vector v);
        
        double distance(Vector v);
        
        float getCoord(int i);
        
        bool toBeErased;
        
    private:
        
    };
}

#endif /* Vector_h */
