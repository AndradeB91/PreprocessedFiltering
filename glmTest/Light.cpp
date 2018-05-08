//
//  Light.cpp
//  glmTest
//
//  Created by Lucas Andrade on 8/4/15.
//  Copyright (c) 2015 Lucas Andrade. All rights reserved.
//

#include "Light.h"

using namespace tdogl;

Light::Light(){}

Light::~Light(){}

void Light::setPosition(glm::vec4 position){
    this->position = position;
}

void Light::setIntensities(glm::vec3 intensities){
    this->intensities = intensities;
}

void Light::setAttenuation(float attenuation){
    this->attenuation = attenuation;
}

void Light::setAmbientCoefficient(float ambientCoefficient){
    this->ambientCoefficient = ambientCoefficient;
}

void Light::setConeAngle(float coneAngle){
    this->coneAngle = coneAngle;
}

void Light::setConeDirection(glm::vec3 coneDirection){
    this->coneDirection = coneDirection;
}

glm::vec4 Light::getPosiion(){
    return this->position;
}
glm::vec3 Light::getIntensities(){
    return this->intensities;
}

float Light::getAttenuation(){
    return this->attenuation;
}

float Light::getAmbientCoefficient(){
    return this->ambientCoefficient;
}

float Light::getConeAngle(){
    return this->coneAngle;
}
glm::vec3 Light::getConeDirection(){
    return this->coneDirection;
}