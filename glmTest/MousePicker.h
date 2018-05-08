//
//  MousePicker.h
//  glmTest
//
//  Created by Lucas Andrade on 14/03/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#ifndef MousePicker_h
#define MousePicker_h

#include <stdio.h>

#include "glm.hpp"
#include "Camera.h"

#include "Mesh.h"


#endif /* MousePicker_h */

namespace Util {
    class MousePicker{
    private:
        glm::vec3 currentRay;
        
        glm::mat4 projectionMatrix;
        glm::mat4 viewMatrix;
        tdogl::Camera camera;
        float width, height;
        Data::Mesh *mesh;
        std::list<int> indexList;
        
        glm::vec3 calculateMouseRay(double mx, double my);
        glm::vec2 getNormalizedDeviceCoords(float mouseX, float mouseY);
        glm::vec4 toEyeCoords(glm::vec4 clipCoords);
        glm::vec3 toWorldCoords(glm::vec4 eyeCoords);
        
    public:
        MousePicker(glm::mat4 projection, float width, float height, Data::Mesh *mesh);
        ~MousePicker();
        glm::vec3 getCurrentRay();
        int update(glm::vec3 camPos, glm::mat4 camViewMatrix, double mx, double my);
        
        bool rayIntersectsTriangle(glm::vec3 p, glm::vec3 d, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2);
        
        int RayTriangleIntersection(glm::vec3 ro, glm::vec3 rd, glm::vec3 va,
                                    glm::vec3 vb, glm::vec3 vc, float *t);
        
        int MeshRayIntersection(glm::vec3 ro, glm::vec3 rd, int &indice);
        
        void addToIndexList(int indice);
        void removeFromIndexList(int indice);
        
        std::list<int> getIndexList();
        
        void indexListClear();
    };
}
