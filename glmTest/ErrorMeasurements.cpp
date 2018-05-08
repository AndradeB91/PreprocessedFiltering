//
//  ErrorMeasurements.cpp
//  glmTest
//
//  Created by Lucas Andrade on 15/03/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#include "ErrorMeasurements.h"

using namespace Util;
using namespace asset;

ErrorMeasurements::ErrorMeasurements(){}

ErrorMeasurements::ErrorMeasurements(asset::ModelAsset groundTruth, asset::ModelAsset filteredModel){
    this->groundTruth = groundTruth;
    this->filteredModel = filteredModel;
}

ErrorMeasurements::ErrorMeasurements(ModelAsset groundTruth, ModelAsset filteredModel,
                                     const std::string &auxFilePath,
                                     const std::string &holeIds){
    this->groundTruth = groundTruth;
    this->filteredModel = filteredModel;
    
    std::ifstream infile(auxFilePath);

    int a;
    bool first = true;
    while(infile >> a){
        if(first){
            this->newVertexAmount = a;
            first = false;
        }
        else this->erasedVertex.push_back(a);
    }
    
    std::ifstream holeFile(holeIds);
    
    int id;
    while(holeFile >> id){
        this->facesIds.push_back(id);
    }
    
    for(int i=0; i<this->facesIds.size(); i++){
        this->faceArray.push_back(groundTruth.getMesh()->faceArray[facesIds[i]]);
    }
}

ErrorMeasurements::~ErrorMeasurements(){}

float ErrorMeasurements::meanDistanceError(){
    
    float sumDist = 0.0;
    int size = static_cast<int>(this->filteredModel.getMesh()->vertexArray.size());
    
    for(int i=0; i<size; i++){
        Vector groundVec   = this->groundTruth.getMesh()->vertexArray[i];
        Vector filteredVec = this->filteredModel.getMesh()->vertexArray[i];
        
        this->initial.push_back(groundVec);
        this->end.push_back(filteredVec);
        
        float res = groundVec.distance(filteredVec);
        sumDist += res;
    }
    
    return sumDist/size;
}

float ErrorMeasurements::meanDistanceToFaces(const std::string &groundTruthPath, bool allMesh){
    
    float sumDist = 0.0;
    int sumDiv = 0;
    int totalSize   = static_cast<int>(this->filteredModel.getMesh()->vertexArray.size());
    int sndPartSize = this->newVertexAmount;
    int fstPartSize = totalSize - sndPartSize;
    
    
    int grThIndex = 0;
    //We use the vertex distance to calculate those vertex that have pairs in the ground-truth...
    for(int i=0; i<fstPartSize; i++){
        
        if(1 == this->erasedVertex[grThIndex]){
            while(1 == this->erasedVertex[grThIndex]){
                grThIndex++;
            }
        }
        
        Vector groundVec   = this->groundTruth.getMesh()->vertexArray[grThIndex];
        Vector filteredVec = this->filteredModel.getMesh()->vertexArray[i];
        
        
        this->initial.push_back(filteredVec);
        this->end.push_back(groundVec);
        float dist = groundVec.distance(filteredVec);
        
        
        sumDist += dist;
        sumDiv++;
        
        grThIndex++;
        
    }
    
    Eigen::MatrixXd groundTruthVertex;
    Eigen::MatrixXi groundTruthFaces;
    
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    
    igl::readOBJ(groundTruthPath, groundTruthVertex, groundTruthFaces);
    
    std::cout << sndPartSize << std::endl;
    //We have all the new mshsurf vertex in this matrix
    Eigen::MatrixXd mshSurfVertex(sndPartSize, 3);
    
    for(int i=fstPartSize; i<totalSize; i++){
        Vector v = this->filteredModel.getMesh()->vertexArray[i];
        
        for(int j=0; j<3; j++){
            mshSurfVertex(i-fstPartSize,j) = (double)v.getCoord(j);
        }
    }
    
    igl::point_mesh_squared_distance(mshSurfVertex, groundTruthVertex, groundTruthFaces, sqrD, I, C);
    
    
    for(int i=0; i<sqrD.size(); i++){
        sumDist += sqrt(sqrD[i]);
        sumDiv++;
    }

    
    if(allMesh){
        for(int i=0; i<sndPartSize; i++){
            Vector v(mshSurfVertex(i,0), mshSurfVertex(i,1), mshSurfVertex(i,2));
            Vector f(C(i,0), C(i,1), C(i,2));
            this->initial.push_back(v);
            this->end.push_back(f);
        }
    }
    
    return sumDist/(float)sumDiv;
}


float ErrorMeasurements::volumeError(){
    float groundTruthVolume = this->groundTruth.getMesh()->volume();
    float filteredMeshVolume = this->filteredModel.getMesh()->volume();
    
    return filteredMeshVolume/groundTruthVolume;
}

glm::vec2 ErrorMeasurements::minMaxFaceQuality(){
    float minQuality = FLT_MAX;
    float maxQuality = -FLT_MAX;
    
    for(int i=0; i<this->filteredModel.getMesh()->faceArray.size(); i++){
        
        Face f = this->filteredModel.getMesh()->faceArray[i];
        float quality = f.quality();
        if(quality < minQuality){
            minQuality = quality;
        }
        if(quality > maxQuality){
            if(quality > 1000){
                f.print();
            }
            maxQuality = quality;
        }
    }
    
    return glm::vec2(minQuality, maxQuality);
}

float ErrorMeasurements::iglError(const std::string &groundTruthPath, const std::string &filteredPath){
    
    /*Eigen::MatrixXd groundTruthVertex, filteredVertex;
    Eigen::MatrixXi groundTruthFaces,  filteredFaces;
    
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    
    igl::readOBJ(groundTruthPath, groundTruthVertex, groundTruthFaces);
    igl::readOBJ(filteredPath,    filteredVertex,    filteredFaces);
    
    igl::point_mesh_squared_distance(filteredVertex, groundTruthVertex, groundTruthFaces, sqrD, I, C);
    
    
    double sum = 0.0;
    for(int i=0; i<sqrD.size(); i++){
        sum += sqrt(sqrD[i]);
    }
    
    //std::cout << "distance: " << sum/sqrD.size() << std::endl;
    
    //for(int i=0; i<filteredVertex.size(); i++){
      //  std::cout << filteredVertex.x() <<  " " << filteredVertex.y() << " " <<  filteredVertex.z() << std::endl;
    //}
    */
    
    Eigen::MatrixXd filteredVertex(2,2);
    
   
    filteredVertex(0,0) = 1;
    filteredVertex(0,1) = 2;
    filteredVertex(1,0) = 3;
    filteredVertex(1,1) = 4;

    
    std::cout << filteredVertex << std::endl;

    
    
    //return sum/sqrD.size();
    return 1;
}









