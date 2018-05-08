//
//  Light.h
//  glmTest
//
//  Created by Lucas Andrade on 8/4/15.
//  Copyright (c) 2015 Lucas Andrade. All rights reserved.
//

#ifndef __glmTest__Light__
#define __glmTest__Light__

#include "glm.hpp"

namespace tdogl {
    class Light{
    public:
        Light();
        ~Light();
        
        void setPosition(glm::vec4 position);
        void setIntensities(glm::vec3 intensities);
        void setAttenuation(float attenuation);
        void setAmbientCoefficient(float ambientCoefficient);
        void setConeAngle(float coneAngle);
        void setConeDirection(glm::vec3 coneDirection);
        
        glm::vec4 getPosiion();
        glm::vec3 getIntensities();
        float getAttenuation();
        float getAmbientCoefficient();
        float getConeAngle();
        glm::vec3 getConeDirection();
        
    private:
        glm::vec4 position;
        glm::vec3 intensities;
        float attenuation;
        float ambientCoefficient;
        float coneAngle;
        glm::vec3 coneDirection;
    };
}



#endif /* defined(__glmTest__Light__) */
