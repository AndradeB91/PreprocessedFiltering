//
//  Ray.cpp
//  glmTest
//
//  Created by Lucas Andrade on 15/03/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#define EPSILON 0.000001

#include "Ray.h"

using namespace Util;
using namespace primitives;

Ray::Ray(glm::vec3 ro, glm::vec3 rd){
    this->ro = ro;
    this->rd = rd;
}

Ray::Ray(){}

Ray::~Ray(){}

void Ray::setRo(glm::vec3 ro){
    this->ro = ro;
}
void Ray::setRd(glm::vec3 rd){
    this->rd = rd;
}

glm::vec3 Ray::getRf(float t){
    return this->ro + t * this->rd;
}

glm::vec3 Ray::getRo(){
    return this->ro;
}

glm::vec3 Ray::getRd(){
    return this->rd;
}
/*
 * Moller-Trumbore intersection.
 * https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
 */
int Ray::RayTriangleIntersection(primitives::Face f, float *t){
    
    glm::vec3 va,vb,vc;

    va.x = f.v0->x; va.y = f.v0->y; va.z = f.v0->z;
    vb.x = f.v1->x; vb.y = f.v1->y; vb.z = f.v1->z;
    vc.x = f.v2->x; vc.y = f.v2->y; vc.z = f.v2->z;
    
    glm::vec3 edge1, edge2, pvec, qvec, tvec;
    edge1 = vb - va;
    edge2 = vc - va;
    
    pvec = glm::cross(this->rd, edge2);
    
    float det = glm::dot(edge1, pvec);
    
    //TODO considerar culling?
    // det < 0 -> backside do triangle
    if (det > -EPSILON && det < EPSILON)
        return 0;
    float invdet = 1.0f/det;
    
    tvec = this->ro - va;
    float u = glm::dot(tvec, pvec)*invdet;
    
    if(u < 0.0f || u > 1.0f)
        return 0;
    
    qvec = glm::cross(tvec, edge1);
    float v = glm::dot(this->rd, qvec)*invdet;
    if(v < 0.0f || u + v > 1.0f)
        return 0;
    
    *t =  glm::dot(edge2, qvec)*invdet;
    if (*t <= EPSILON)
        return 0;
    
    return 1;
}

bool Ray::RayTriangleIntersection2(primitives::Face face, float *t) {
    
    glm::vec3 v0,v1,v2;
    
    v0.x = face.v0->x; v0.y = face.v0->y; v0.z = face.v0->z;
    v1.x = face.v1->x; v1.y = face.v1->y; v1.z = face.v1->z;
    v2.x = face.v2->x; v2.y = face.v2->y; v2.z = face.v2->z;
    
    glm::vec3 e1,e2,h,s,q;
    float a,f,u,v;
    //float t;
    
    e1 = v1 - v0;
    e2 = v2 - v0;
    
    h = glm::cross(this->rd, e2);
    
    a = glm::dot(e1,h);
    
    if (a > -0.00001 && a < 0.00001)
        return false;
    
    f = 1/a;
    s = this->ro - v0;
    u = f * (glm::dot(s,h));
    
    if (u < 0.0 || u > 1.0)
        return false;
    
    q = glm::cross(s, e1);
    v = f * glm::dot(this->rd, q);
    
    if (v < 0.0 || u + v > 1.0)
        return false;
    
    // at this stage we can compute t to find out where
    // the intersection point is on the line
    *t = f * glm::dot(e2,q);
    
    if (*t > 0.00001) // ray intersection
        return true;
    
    else // this means that there is a line intersection
        // but not a ray intersection
        return false;
    
}

bool Ray::RayFacesIntersection(std::vector<primitives::Face> faceArray, float *dist){
    bool intersect = false;
    float minT = FLT_MAX;
    int faceTouched = -1;
    
    //Run for each face in mesh
    for(int i=0; i<faceArray.size(); i++){
        
        Face f = faceArray[i];
        float newt;
        
        Vector n = f.normal.normalize();
        this->rd.x = n.x;
        this->rd.y = n.y;
        this->rd.z = n.z;
        
        if(this->RayTriangleIntersection2(f, &newt)){
            if(newt < minT){
                minT = newt;
                intersect = true;
                faceTouched = i;
            }
        }
        
        this->rd.x = -n.x;
        this->rd.y = -n.y;
        this->rd.z = -n.z;
        
        if(this->RayTriangleIntersection2(f, &newt)){
            if(newt < minT){
                minT = newt;
                intersect = true;
                faceTouched = i;
            }
        }
    }
    
    //Now we calculate the distance between ro and intersection point
    if(true == intersect){
        glm::vec3 rf = this->getRf(minT);
        *dist = glm::distance(this->ro, rf);
        //std::cout << *dist << " " << faceTouched << std::endl;
        //std::cout << faceTouched << ",";
    }
    return intersect;
}


















