//
//  ModelAsset.h
//  glmTest
//
//  Created by Lucas Andrade on 7/14/15.
//  Copyright (c) 2015 Lucas Andrade. All rights reserved.
//

#ifndef __glmTest__ModelAsset__
#define __glmTest__ModelAsset__

#include <stdio.h>
#include "Program.h"
#include "Texture.h"
#include "Bitmap.h"

#include "Face.h"
#include "Vector.h"
#include "Mesh.h"
#include "MainDrive.h"

#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing fla

#include "glm/glm/gtc/matrix_transform.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <math.h>       /* exp */



namespace asset {
    class ModelAsset{
    public:
        ModelAsset();
             
        
        ModelAsset(aiMesh *mesh,
                   const std::string &vertexShader,
                   const std::string &fragmentShader,
                   const std::string &texturePath,
                   const GLenum drawType,
                   const bool light,
                   glm::mat4 camera,
                   glm::mat4 &transformation,
                   GLfloat shininess,
                   glm::vec3 specularColor,
                   bool hasUV);
        
        ModelAsset(const std::string &objFilePath,
                   const std::string &vertexShader,
                   const std::string &fragmentShader,
                   const std::string &texturePath,
                   const std::string &edgeTexturePath,
                   const GLenum drawType,
                   const bool light,
                   glm::mat4 camera,
                   glm::mat4 &transformation,
                   GLfloat shininess,
                   glm::vec3 specularColor);
        
        ModelAsset(Data::Mesh *m,
                   const std::string &vertexShader,
                   const std::string &fragmentShader,
                   const std::string &texturePath,
                   const std::string &edgeTexturePath,
                   const GLenum drawType,
                   const bool light,
                   glm::mat4 camera,
                   glm::mat4 &transformation,
                   GLfloat shininess,
                   glm::vec3 specularColor);

        
        ~ModelAsset();
        
        void readBaseMesh(const std::string &objFilePath);
        
        tdogl::Program *getShaders();
        tdogl::Texture *getTexture();
        tdogl::Texture *getEdgeTexture();

        
        glm::mat4 *getCameraMatrix();

        glm::mat4 *getTransformationMatrix();
        
        bool getLightEnable();
        
        
        //GETS...
        GLuint getVbo();
        GLuint getVao();
        GLuint getNormalVbo();
        GLuint getNormalVao();
        GLuint getEdgeVbo();
        GLuint getEdgeVao();
        GLuint getDrawType();
        GLuint getDrawStart();
        GLuint getDrawCount();
        GLuint getDrawEdgeCount();

        
        GLfloat getShininess();
        glm::vec3 getSpecularColor();
        unsigned long getBufferSize();
        
        tdogl::Program *LoadShaders(const char *vertFileName,
                                    const char *fragFileName);
        
        void LoadTexture(const std::string &texturePath);
        
        void LoadProjectionMatrix(glm::mat4 &projection);
        void LoadCameraMatrix(glm::mat4 &camera);
        
        
    
        
        void makeEdgeConnections(Face *f);
        
        //methods for the consistency function
        double phi(std::list<Face *> P);
        double R(std::list<Face *> P);
        double H(std::list<Face *> P);
        std::list<Face *> edgeNeighbor(std::list<Face *> P, Face &f);
        Vector patchAverageNormal(std::list<Face *> P);
        
        
        void computeMeshGuidances();
        void updatingVertices();
        void updateNormals();
        void updateEdges();
        void updatingNormalsFromTriangles();
        
        void showMshSurfMesh();
        void showFinalMesh();
        void showMeshWithHoles();
        void showOriginalMesh();
        
        Vector filteredNormal(Face &f);
        //kernels
        double spatialKernel(Vector p, Vector q);
        double rangeKernel(Vector Ip, Vector Iq);

        void processAllMesh();
        void executeMshSurf();
        void executeMshSurf(std::list<int> indexList);
        Data::Mesh *getMesh();
        Data::Mesh *getMshSurfMesh();
        Data::Mesh *getOriginalMesh();
        int getNewVertexAmount();
        
        void putInBuffer(GLuint &vbo, GLuint &vao, GLfloat data[], unsigned long size);

        void renderSelectedFaces(std::list<int> indexList);
        void renderSelectedVertex(std::vector<Vector> initial, std::vector<Vector> end);
        void renderMesh(Data::Mesh *mesh);
        int joinMeshesData(Data::Mesh *newMesh);
        
        int findIndex(Face &f);
        //bool isInPatch();
        std::list<Face *> patchToBeDeleted();
        std::list<Face *> patchToBeDeleted(std::list<int> indexList);
        Data::Mesh openHole();
        Data::Mesh openHole(std::list<int> indexList);
        std::vector<Face> removeToBeErasedFaces();
        void removeToBeErasedVertex();
        void removeToBeErasedEdges();
        void removeToBeErasedNormals();

    private:
        
        MeshSurf::MainDrive mainDrive;

        Data::Mesh originalMesh;
        Data::Mesh mesh; //The processed final mesh.
        Data::Mesh mshSurf;
        Data::Mesh meshWithHoles;
        Data::Mesh baseMeshForMshSurf;
        
        int newVertexAmount;

        tdogl::Program *shaders;
        tdogl::Texture *texture;
        tdogl::Texture *edgeTexture;

        
        std::vector<tdogl::Texture*> textures;
        
        glm::mat4 *camera;        
        glm::mat4 *transformation;
        
        bool lightEnable;
        
        GLuint vbo;
        GLuint vao;
        
        GLuint normalVbo;
        GLuint normalVao;
        
        GLuint edgeVbo;
        GLuint edgeVao;
        
        GLenum drawType;
        GLint drawStart;
        GLint drawCount;
        GLint drawEdgeCount;
        
        GLfloat shininess;
        glm::vec3 specularColor;
        
        
    };
}



#endif /* defined(__glmTest__ModelAsset__) */
