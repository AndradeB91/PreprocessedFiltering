//
//  MousePicker.cpp
//  glmTest
//
//  Created by Lucas Andrade on 14/03/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#define EPSILON 0.000001

#include "MousePicker.h"

using namespace Util;

MousePicker::MousePicker(glm::mat4 projection, float width, float height, Data::Mesh *mesh){
    this->projectionMatrix = projection;
    this->width = width;
    this->height = height;
    this->mesh = mesh;
}

MousePicker::~MousePicker(){};

glm::vec3 MousePicker::calculateMouseRay(double mx, double my){
    float mouseX = mx;
    float mouseY = my;
    glm::vec2 normalizedCoords = this->getNormalizedDeviceCoords(mouseX, mouseY);
    glm::vec4 clipCoords = glm::vec4(normalizedCoords.x, normalizedCoords.y, -1.0f, 1.0f);
    glm::vec4 eyeCoords = this->toEyeCoords(clipCoords);
    glm::vec3 worldRay = this->toWorldCoords(eyeCoords);
    return worldRay;
}

glm::vec2 MousePicker::getNormalizedDeviceCoords(float mouseX, float mouseY){
    float x = (2.0f * mouseX) / this->width - 1.0f;
    float y = (2.0f * mouseY) / this->height - 1.0f;
    glm::vec2 res(x,-y);
    return res;
}

glm::vec4 MousePicker::toEyeCoords(glm::vec4 clipCoords){
    glm::mat4 invertedProjection = glm::inverse(this->projectionMatrix);
    glm::vec4 eyeCoords = invertedProjection * clipCoords;
    glm::vec4 res(eyeCoords.x, eyeCoords.y, -1.0f, 0.0f);
    return res;
}

glm::vec3 MousePicker::toWorldCoords(glm::vec4 eyeCoords){
    glm::mat4 invertedView = glm::inverse(this->viewMatrix);
    glm::vec4 rayWorld = invertedView * eyeCoords;
    glm::vec3 mouseRay = glm::vec3(rayWorld.x, rayWorld.y, rayWorld.z);
    return glm::normalize(mouseRay);
}

int MousePicker::update(glm::vec3 camPos, glm::mat4 camViewMatrix, double mx, double my){
    this->viewMatrix = camViewMatrix;
    this->currentRay = calculateMouseRay(mx, my);
    
    int indice = -1;
    this->MeshRayIntersection(camPos, this->currentRay, indice);
    return indice;
}

int MousePicker::MeshRayIntersection(glm::vec3 ro, glm::vec3 rd, int &indice){
    
    int intersect = 0;
    float minT = FLT_MAX;
    
    //Run for each face in mesh
    for(int i=0; i<mesh->faceArray.size(); i++){
        Face *f = &mesh->faceArray[i];
        glm::vec3 va, vb, vc;
        va.x = f->v0->x; va.y = f->v0->y; va.z = f->v0->z;
        vb.x = f->v1->x; vb.y = f->v1->y; vb.z = f->v1->z;
        vc.x = f->v2->x; vc.y = f->v2->y; vc.z = f->v2->z;
        
        
        float newt;
        if(RayTriangleIntersection(ro, rd, va, vb, vc, &newt)){
            if(newt < minT){
                minT = newt;
                intersect = 1;
                indice = f->id;
            }
        }
    }
    return intersect;
}

glm::vec3 MousePicker::getCurrentRay(){
    return this->currentRay;
}

void MousePicker::addToIndexList(int indice){
    bool in = false;
    for(std::list<int>::iterator i = this->indexList.begin(); i != this->indexList.end(); i++){
        if(*i == indice) in = true;
    }
    if(false == in && -1 != indice)
        this->indexList.push_back(indice);
}

void MousePicker::removeFromIndexList(int indice){
    std::list<int>::iterator i;
    for(i = this->indexList.begin(); i != this->indexList.end(); i++){
        if(*i == indice){
            this->indexList.erase(i);
        }
    }
}

std::list<int> MousePicker::getIndexList(){
    return this->indexList;
}

void MousePicker::indexListClear(){
    this->indexList.clear();
}

bool MousePicker::rayIntersectsTriangle(glm::vec3 p, glm::vec3 d, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2) {
    
    glm::vec3 e1,e2,h,s,q;
    float a,f,u,v;
    float t;
    
    e1 = v1 - v0;
    e2 = v2 - v0;
    
    h = glm::cross(d, e2);

    a = glm::dot(e1,h);
    
    if (a > -0.00001 && a < 0.00001)
        return false;
    
    f = 1/a;
    s = p - v0;
    u = f * (glm::dot(s,h));
    
    if (u < 0.0 || u > 1.0)
        return false;
    
    q = glm::cross(s, e1);
    v = f * glm::dot(d, q);
    
    if (v < 0.0 || u + v > 1.0)
        return false;
    
    // at this stage we can compute t to find out where
    // the intersection point is on the line
    t = f * glm::dot(e2,q);
    
    if (t > 0.00001) // ray intersection
        return true;
    
    else // this means that there is a line intersection
        // but not a ray intersection
        return false;
    
}

/*
 * Moller-Trumbore intersection.
 * https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
 */
int MousePicker::RayTriangleIntersection(glm::vec3 ro, glm::vec3 rd, glm::vec3 va,
                                     glm::vec3 vb, glm::vec3 vc, float *t){
    
    glm::vec3 edge1, edge2, pvec, qvec, tvec;
    edge1 = vb - va;
    edge2 = vc - va;
    
    pvec = glm::cross(rd, edge2);
    
    float det = glm::dot(edge1, pvec);
    
    //TODO considerar culling?
    // det < 0 -> backside do triangle
    if (det > -EPSILON && det < EPSILON)
        return 0;
    float invdet = 1.0f/det;
    
    tvec = ro - va;
    float u = glm::dot(tvec, pvec)*invdet;
    
    if(u < 0.0f || u > 1.0f)
        return 0;
    
    qvec = glm::cross(tvec, edge1);
    float v = glm::dot(rd, qvec)*invdet;
    if(v < 0.0f || u + v > 1.0f)
        return 0;
    
    *t =  glm::dot(edge2, qvec)*invdet;
    if (*t <= EPSILON)
        return 0;
    
    return 1;
}


