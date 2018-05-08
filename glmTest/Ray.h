//
//  Ray.h
//  glmTest
//
//  Created by Lucas Andrade on 15/03/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#ifndef Ray_h
#define Ray_h

#include <stdio.h>
#include <vector>

#include "glm.hpp"
#include "Face.h"


namespace Util {
    class Ray{
    private:
        glm::vec3 ro;
        glm::vec3 rd;
    public:
        Ray(glm::vec3 ro, glm::vec3 rd);
        Ray();
        ~Ray();
        void setRo(glm::vec3 ro);
        void setRd(glm::vec3 rd);
        glm::vec3 getRo();
        glm::vec3 getRd();
        glm::vec3 getRf(float t);
        bool RayTriangleIntersection2(primitives::Face face, float *t);
        int RayTriangleIntersection(primitives::Face f, float *t);
        bool RayFacesIntersection(std::vector<primitives::Face> faceArray, float *dist);
    };
}

#endif /* Ray_h */
