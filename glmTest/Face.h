//
//  Face.hpp
//  glmTest
//
//  Created by Lucas Andrade on 30/09/16.
//  Copyright © 2016 Lucas Andrade. All rights reserved.
//


#ifndef Face_hpp
#define Face_hpp

#include <stdio.h>

#include <assimp/scene.h>           // Output data structure
#include "glm/glm/gtc/matrix_transform.hpp"

#include "Vector.h"
#include "Edge.h"

#include <iostream>

namespace primitives {
        
    class Face{
    public:
        Face();
        ~Face();
        
        int id;
        bool processed;
        Vector *v0, *v1, *v2;
        Edge *e0, *e1, *e2;
        //glm::vec3 normal;
        Vector normal;
        
        double consistency;
        Vector guidance;
        bool toBeErased;
        
        float area();
        Vector center();
        Vector *getVector(int i);
        
        std::list<Face *> getFacePatch();
        std::list<Face *> realFaceNeighbors();
        
        float signedVolume();
        float quality();
        void print();
    private:
        
    };
}


#endif /* Face_h */
