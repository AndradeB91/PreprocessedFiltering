//
//  MainDrive.hpp
//  glmTest
//
//  Created by Lucas Andrade on 26/01/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#ifndef MainDrive_h
#define MainDrive_h

#include <stdio.h>
#include <unordered_map>
#include <iostream>

//This file has the functions: "MshSurfSetSupportMesh" and "MshSurf3D"
#include "tecgraf/Surf3D/mshsurf3d.h"

//#include "ModelAsset.h"
#include "Mesh.h"

namespace MeshSurf {
    
    class MainDrive{
        
    public:
        MainDrive();
        ~MainDrive();
        
        void setBaseMesh(Data::Mesh &baseMesh);
        void setMesh(Data::Mesh &mesh);
        void setSupportComputeCurvature(bool b);
        void setSupportMaxElementSize(double maxElemSize);
        void setSupportMaxElementSize();
        void setInvertEdges(bool invert);
        
        void ExecuteSupport();
        
        Data::Mesh getGeneratedMesh();

    private:
        Data::Mesh *mesh;
        Data::Mesh generatedMesh;
        
        Data::Mesh *baseMesh;
        
        bool invert;
        
        Surf3DMessFunc *messageFunction;
        bool supportComputeCurvature;
        double supportMaxElementSize;
        
        void convertToMeshSurf(int &nPts, double *&pts, int &nElems, int *&elems,
                               Data::Mesh *mesh, int type, bool clockWise);
        
        bool convertToBoundarySurf(int &nPts, double *&pts, int &nBoundEdges, int &nInterEdges, int *&edges,
                                  std::unordered_map<int, Vector> &map, Data::Boundary *boundary, bool clockWise);
        
        bool convertToMesh(int nPts, double *pts, int nElems, int *elems, int type,
                           std::unordered_map<int, Vector> &map, bool clockWise, bool rebuildMap, bool checkBoundary);
        
    };
}

#endif /* MainDrive_h */
