//
//  MainDrive.cpp
//  glmTest
//
//  Created by Lucas Andrade on 26/01/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#include "MainDrive.h"



using namespace MeshSurf;


MainDrive::MainDrive(){
    this->messageFunction = NULL;
    this->supportMaxElementSize = 0.0;
    this->supportComputeCurvature = false;
}

MainDrive::~MainDrive(){};

void MainDrive::setBaseMesh(Data::Mesh &baseMesh){
    this->baseMesh = &baseMesh;
}

void MainDrive::setMesh(Data::Mesh &mesh){
    this->mesh = &mesh;
}

void MainDrive::setSupportComputeCurvature(bool b){
    this->supportComputeCurvature = b;
}

void MainDrive::setSupportMaxElementSize(double maxElemSize){
    this->supportMaxElementSize = maxElemSize;
}

void MainDrive::setSupportMaxElementSize(){
    
    double meanDistance = 0.0;
    std::list<Edge>::iterator boundEdge;
    
    for(boundEdge = this->mesh->boundary.outline.begin();
        boundEdge != this->mesh->boundary.outline.end(); boundEdge++){
        
        Vector *v0 = boundEdge->v0;
        Vector *v1 = boundEdge->v1;
        
        meanDistance = v0->distance(*v1);
    }
    
    meanDistance /= this->mesh->boundary.outline.size();
    
    //The mean edge distance is too small so we increase it by a factor
    this->supportMaxElementSize = meanDistance * 2.5;
}

void MainDrive::setInvertEdges(bool invert){
    this->invert = invert;
}

void MainDrive::ExecuteSupport(){
    
    
    int    nSupPts = 0;
    double *supPts = NULL;
    int    nSupElems = 0;
    int    *supElems = NULL;
    
    //We have to pass the base mesh
    if(this->baseMesh != NULL){
        this->convertToMeshSurf(nSupPts, supPts, nSupElems, supElems, this->baseMesh, 3, false);
    }else{
        this->convertToMeshSurf(nSupPts, supPts, nSupElems, supElems, this->mesh, 3, false);
    }
    
    
    int    nPts = 0;
    double *pts = NULL;
    int    nBoundEdges = 0;
    int    nInterEdges = 0;
    int    *edges = NULL;
 
    std::unordered_map<int, Vector> map;
    
    //Now we pass the boundary for the mshsurf to be created

    if(this->invert){
        for(std::list<Edge>::iterator it = this->mesh->boundary.outline.begin(); it != this->mesh->boundary.outline.end(); it++){

            Vector *aux = it->v0;
            it->v0 = it->v1;
            it->v1 = aux;
        }
    }


    this->convertToBoundarySurf(nPts, pts, nBoundEdges, nInterEdges, edges, map, &this->mesh->boundary, false);

    
    //Right below we are using rmsh functions: MshSurfSetSupportMesh() and MshSurf3D()
    
    MshSurfSetSupportMesh(NULL, nSupPts, supPts, nSupElems, supElems);
    

    int nGenPts = 0;
    double *genPts = NULL;
    int nGenElems = 0;
    int *genElems = NULL;
    
    MshSurf3D(nPts, pts, nBoundEdges, nInterEdges, edges,
              this->supportMaxElementSize, this->supportComputeCurvature ? 1 : 0, this->messageFunction,
              &nGenPts, &genPts, &nGenElems, &genElems);
    
    
    this->convertToMesh(nGenPts, genPts, nGenElems, genElems, 0, map, false, false, false);
    
    delete [] supPts;
    delete [] supElems;
    delete [] pts;
    delete [] edges;
    free(genPts);
    free(genElems);
}

void MainDrive::convertToMeshSurf(int &nPts, double *&pts, int &nElems, int *&elems,
                                  Data::Mesh *mesh, int type, bool clockWise){
    
    
    nPts = static_cast<int>(mesh->vertexArray.size());
    nElems = static_cast<int>(mesh->faceArray.size());
    
    pts = new double[3*nPts];
    elems = new int[(type+1)*nElems];
    
    std::unordered_map<int,int> vMap;
    
    for(int i=0; i<nPts; i++){
        pts[3*i+0] = mesh->vertexArray[i].x;
        pts[3*i+1] = mesh->vertexArray[i].y;
        pts[3*i+2] = mesh->vertexArray[i].z;
        
        vMap[mesh->vertexArray[i].id] = i;
    }
    
    for(int i=0; i<nElems; i++){
        elems[(type+1)*i + 0] = type;
        
        Face f = mesh->faceArray[i];
        for(int j=0; j<type; j++){
            
            int jj = clockWise ? (type - j) % type : j;
            elems[(type+1)*i + j+1] = vMap[f.getVector(jj)->id];
        }
        
    }
}

bool MainDrive::convertToBoundarySurf(int &nPts, double *&pts, int &nBoundEdges, int &nInterEdges, int *&edges,
                                      std::unordered_map<int, Vector> &map,  Data::Boundary *boundary,
                                      bool clockWise){
    
    nPts = static_cast<int>(boundary->outline.size());
    nBoundEdges = static_cast<int>(boundary->outline.size());
    nInterEdges = 0;
    
    pts = new double[3*nPts];
    edges = new int[2*nBoundEdges];
    
    std::unordered_map<int, int> vMap;
    
    map.reserve(nPts);

    std::list<Edge>::iterator bound = boundary->outline.begin();
    for(int i=0; i<nPts; i++){
        
        Vector *v0 = bound->v0;
        
        pts[3*i+0] = v0->x;
        pts[3*i+1] = v0->y;
        pts[3*i+2] = v0->z;
        
        vMap[v0->id] = i;
        
        bound++;
    }
    std::list<Edge>::iterator eiter = boundary->outline.begin();
    int j = clockWise ? 1 : 0;
    
    for (int i=0; i<nBoundEdges; i++){
        
        Vector *v1 = eiter->getVector(j);
        Vector *v2 = eiter->getVector(1-j);
        
        edges[2*i+0] = vMap[v1->id];
        edges[2*i+1] = vMap[v2->id];
        eiter++;
    }
    
    return true;
}

bool MainDrive::convertToMesh(int nPts, double *pts, int nElems, int *elems, int type,
                              std::unordered_map<int, Vector> &map, bool clockWise, bool rebuildMap, bool checkBoundary){
    
    if ((nPts == 0) || (nElems == 0) || (pts == NULL) || (elems == NULL))
    {
        return false;
    }

    this->generatedMesh.convertFromMshSurfToMesh(nPts, pts, nElems, elems);
    return true;
}

Data::Mesh MainDrive::getGeneratedMesh(){
    return this->generatedMesh;
}




















