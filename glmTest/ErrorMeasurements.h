//
//  ErrorMeasurements.h
//  glmTest
//
//  Created by Lucas Andrade on 15/03/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#ifndef ErrorMeasurements_h
#define ErrorMeasurements_h

#include <stdio.h>

#include "ModelAsset.h"
#include "Ray.h"

#include <libigl/include/igl/point_mesh_squared_distance.h>
#include <libigl/include/igl/readOBJ.h>


namespace Util {
    class ErrorMeasurements{
    private:
        asset::ModelAsset groundTruth;
        asset::ModelAsset filteredModel;
        
        int newVertexAmount;
        std::vector<int> erasedVertex;
        
        std::vector<int> facesIds;
        std::vector<Face> faceArray;
        
    public:
        
        std::vector<Vector> initial, end;
        
        ErrorMeasurements();
        ErrorMeasurements(asset::ModelAsset groundTruth, asset::ModelAsset filteredModel);
        ErrorMeasurements(asset::ModelAsset groundTruth, asset::ModelAsset filteredModel,
                          const std::string &auxFilePath,
                          const std::string &holeIds);
        ~ErrorMeasurements();
        
        float meanDistanceError();
        
        float meanDistanceToFaces(const std::string &groundTruthPath, bool allMesh);
        
        float volumeError();
        
        glm::vec2 minMaxFaceQuality();
        
        float iglError(const std::string &groundTruthPath, const std::string &filteredPath);
    };
}

#endif /* ErrorMeasurements_h */
