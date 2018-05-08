//
//  Scene.h
//  glmTest
//
//  Created by Lucas Andrade on 1/18/16.
//  Copyright (c) 2016 Lucas Andrade. All rights reserved.
//

#ifndef __glmTest__Scene__
#define __glmTest__Scene__

#include <stdio.h>
#include <vector>

#include "Program.h"
#include "Texture.h"
#include "Bitmap.h"
#include "ModelAsset.h"

#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla

#include "glm/glm/gtc/matrix_transform.hpp"

namespace asset {
    class Scene{
    public:
        Scene(const std::string &objFilePath,
              const std::string &vertexShader,
              const std::string &fragmentShader,
              const GLenum drawType,
              const bool light,
              glm::mat4 camera,
              glm::mat4 &transformation);
        
        std::vector<ModelAsset> getAssets();
    
    private:
        std::vector<ModelAsset> assets;
        
    };
}

#endif /* defined(__glmTest__Scene__) */
