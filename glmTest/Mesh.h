//
//  Mesh.hpp
//  glmTest
//
//  Created by Lucas Andrade on 30/01/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#ifndef Mesh_h
#define Mesh_h

#include <stdio.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <unordered_map>

#include "Boundary.h"

#include "Face.h"
#include "Vector.h"
#include "Edge.h"

using namespace primitives;

namespace Data{
    
    class Mesh{
    public:
        Mesh();
        ~Mesh();
        
        void setName(std::string name);
        std::string getName();
        
        std::vector<Face>   faceArray;
        std::vector<Vector> vertexArray;
        
        std::list<Vector> normalList;
        std::list<Edge> edgeList;
        
        Boundary boundary;
        
        float volume();

        void findOutline();
        

        void makeVertexArrayFromFaceArray();
        void makeVertexArrayFromDuplicatedFaceArray();
        int countReps(Vector v);
        
        void remakeFacePointers();

        
        bool isInVertexArray(Vector v);
        
        bool constructDataStructureForHoleMesh();
        bool convertFromAssimpToMesh(aiMesh *assimpMesh, std::string objFilePath);
        bool convertFromAssimpToMesh_simpleVersion(aiMesh *assimpMesh);
        bool convertFromMshSurfToMesh(int nPts, double *pts, int nElems, int *elems);
        
        int at(aiVector3D vec);
        int at(Vector vec);
        
        void makeEdgeConnections(Face *f);
        
        bool getHasNormals();
        bool getHasTextureCoords();
        
        bool exportObj(int newVertexAmount, std::string *name);
        bool exportHoleFacesId(std::list<int> facesArrayIndex);
        
    private:
        bool hasNormals;
        bool hasTextureCoords;
        
        std::string name;
    };
    
}

#endif /* Mesh_h */
