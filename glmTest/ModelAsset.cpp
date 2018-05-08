//
//  ModelAsset.cpp
//  glmTest
//
//  Created by Lucas Andrade on 7/14/15.
//  Copyright (c) 2015 Lucas Andrade. All rights reserved.
//

#include "ModelAsset.h"

#define NORMAL_FACTOR 1.2;


using namespace asset;
using namespace primitives;
using namespace std;
using namespace Data;
//using namespace MeshSurf;

ModelAsset::ModelAsset(){}


ModelAsset::ModelAsset(aiMesh *mesh,
                       const std::string &vertexShader,
                       const std::string &fragmentShader,
                       const std::string &texturePath,
                       const GLenum drawType,
                       const bool light,
                       glm::mat4 camera,
                       glm::mat4 &transformation,
                       GLfloat shininess,
                       glm::vec3 specularColor,
                       bool hasUV){
    
    
    int numVerts = mesh->mNumFaces*3;
    GLfloat vertexData[ numVerts*8 ];
    int k = 0;
    
    for(unsigned int i=0; i<mesh->mNumFaces; i++){
        const aiFace& face = mesh->mFaces[i];
        
        for(int j=0; j<3; j++){
            
            aiVector3D pos = mesh->mVertices[face.mIndices[j]];
            
            vertexData[k++] = pos.x;
            vertexData[k++] = pos.y;
            vertexData[k++] = pos.z;
            
            if(!mesh->HasTextureCoords(0)){
                vertexData[k++] = 0.0;
                vertexData[k++] = 0.0;
            }else{
                aiVector3D uv = mesh->mTextureCoords[0][face.mIndices[j]];
                vertexData[k++] = uv.x;
                vertexData[k++] = uv.y;
            }
            
            
            aiVector3D normal = mesh->mNormals[face.mIndices[j]];
            vertexData[k++] = normal.x;
            vertexData[k++] = normal.y;
            vertexData[k++] = normal.z;
            
        }
    }
    
    const char *vshader = strcpy((char*)malloc(vertexShader.length()+1), vertexShader.c_str());
    const char *fshader = strcpy((char*)malloc(fragmentShader.length()+1), fragmentShader.c_str());
    
    this->shaders = LoadShaders(vshader, fshader);
    this->drawType = drawType;
    this->drawStart = 0;
    this->drawCount = numVerts;
    
    LoadTexture(texturePath);
    LoadCameraMatrix(camera);
    this->transformation = &transformation;
    if(light) this->lightEnable = true;
    this->shininess = shininess;
    this->specularColor = specularColor;
    
    this->putInBuffer(this->vbo, this->vao, vertexData, sizeof(vertexData));

}


ModelAsset::ModelAsset(Data::Mesh *m,
                       const std::string &vertexShader,
                       const std::string &fragmentShader,
                       const std::string &texturePath,
                       const std::string &edgeTexturePath,
                       const GLenum drawType,
                       const bool light,
                       glm::mat4 camera,
                       glm::mat4 &transformation,
                       GLfloat shininess,
                       glm::vec3 specularColor){
    
    const char *vshader = strcpy((char*)malloc(vertexShader.length()+1), vertexShader.c_str());
    const char *fshader = strcpy((char*)malloc(fragmentShader.length()+1), fragmentShader.c_str());
    
    this->mesh = *m;
    
    this->shaders = LoadShaders(vshader, fshader);
    this->drawType = drawType;
    this->drawStart = 0;
    
    LoadTexture(texturePath);
    LoadCameraMatrix(camera);
    this->transformation = &transformation;
    if(light) this->lightEnable = true;
    this->shininess = shininess;
    this->specularColor = specularColor;
    
    tdogl::Bitmap bitmap = tdogl::Bitmap::bitmapFromFile(edgeTexturePath);
    bitmap.flipVertically();
    this->edgeTexture = new tdogl::Texture(bitmap);
    textures.push_back(new tdogl::Texture(bitmap));
    
    this->renderMesh(&this->mesh);
}

ModelAsset::ModelAsset(const std::string &objFilePath,
                       const std::string &vertexShader,
                       const std::string &fragmentShader,
                       const std::string &texturePath,
                       const std::string &edgeTexturePath,
                       const GLenum drawType,
                       const bool light,
                       glm::mat4 camera,
                       glm::mat4 &transformation,
                       GLfloat shininess,
                       glm::vec3 specularColor){

    this->mesh.setName(objFilePath);
    
    const char *vshader = strcpy((char*)malloc(vertexShader.length()+1), vertexShader.c_str());
    const char *fshader = strcpy((char*)malloc(fragmentShader.length()+1), fragmentShader.c_str());
    
    Assimp::Importer importer;
    const char *file_path = strcpy((char*)malloc(objFilePath.length()+1), objFilePath.c_str());
    const aiScene *scene = importer.ReadFile(file_path, aiProcess_JoinIdenticalVertices);
    aiMesh *assimpMesh = scene->mMeshes[0];
    
    this->mesh.convertFromAssimpToMesh(assimpMesh, objFilePath);
    
//    for(int i=0; i<this->mesh.vertexArray.size();i++){
//        Vector v = this->mesh.vertexArray[i];
//        std::cout << v.x << " " << v.y << " " << v.z << std::endl;
//    }
//    std::cout << "\n" << std::endl;
    //this->mesh.convertFromAssimpToMesh_simpleVersion(assimpMesh);
    this->originalMesh = mesh;
    

    this->shaders = LoadShaders(vshader, fshader);
    this->drawType = drawType;
    this->drawStart = 0;

    LoadTexture(texturePath);
    LoadCameraMatrix(camera);
    this->transformation = &transformation;
    if(light) this->lightEnable = true;
    this->shininess = shininess;
    this->specularColor = specularColor;
    
    tdogl::Bitmap bitmap = tdogl::Bitmap::bitmapFromFile(edgeTexturePath);
    bitmap.flipVertically();
    this->edgeTexture = new tdogl::Texture(bitmap);
    textures.push_back(new tdogl::Texture(bitmap));
    
    //this->processAllMesh();
        
    this->renderMesh(&this->originalMesh);
}

ModelAsset::~ModelAsset(){}

void ModelAsset::readBaseMesh(const std::string &objFilePath){
    Assimp::Importer importer;
    const char *file_path = strcpy((char*)malloc(objFilePath.length()+1), objFilePath.c_str());
    const aiScene *scene = importer.ReadFile(file_path, aiProcess_JoinIdenticalVertices);
    aiMesh *assimpMesh = scene->mMeshes[0];
    this->baseMeshForMshSurf.convertFromAssimpToMesh(assimpMesh, objFilePath);
}

void ModelAsset::showMshSurfMesh(){
    this->renderMesh(&this->mshSurf);
}

void ModelAsset::showFinalMesh(){
    this->renderMesh(&this->mesh);
}

void ModelAsset::showMeshWithHoles(){
    this->renderMesh(&this->meshWithHoles);
}

void ModelAsset::showOriginalMesh(){
    this->renderMesh(&this->originalMesh);
}


void ModelAsset::processAllMesh(){
    int processIterations = 1;
    for(int i=0; i<processIterations; i++){
        this->executeMshSurf();
    }
}

void ModelAsset::executeMshSurf(std::list<int> indexList){
    Data::Mesh hole = this->openHole(indexList);
    Data::Mesh *mesh = &hole;
    mesh->constructDataStructureForHoleMesh();
    this->mainDrive.setMesh(*mesh);
    this->mainDrive.setInvertEdges(false);
    this->mainDrive.setSupportComputeCurvature(true);
    this->mainDrive.setSupportMaxElementSize(0.07);
    this->mainDrive.ExecuteSupport();
    this->mshSurf = mainDrive.getGeneratedMesh();
    
    std::cout << this->mshSurf.faceArray.size() << std::endl;

    //this->mesh = mshSurf;
    this->meshWithHoles = this->mesh;
    this->newVertexAmount = this->joinMeshesData(&this->mshSurf);
    this->mesh.exportObj(this->newVertexAmount, NULL);
    
//    this->mesh.exportObj(-1, NULL);
//    std::string meshName = "hole_";
//    this->meshWithHoles.exportObj(0, &meshName);
//    this->mesh.exportHoleFacesId(indexList);
}

void ModelAsset::executeMshSurf(){
    
    this->mainDrive.setMesh(this->mesh);
    this->mainDrive.setBaseMesh(this->baseMeshForMshSurf);
    this->mainDrive.setInvertEdges(true);
    this->mainDrive.setSupportComputeCurvature(true);
    this->mainDrive.setSupportMaxElementSize(0.05);
    this->mainDrive.ExecuteSupport();
    this->mshSurf = mainDrive.getGeneratedMesh();
    
    this->meshWithHoles = this->mesh;
    this->newVertexAmount = this->joinMeshesData(&this->mshSurf);
    
    this->mesh.exportObj(this->newVertexAmount, NULL);
}

int ModelAsset::findIndex(Face &f){
    int index = -1;
    for(int i=0; i<this->mesh.faceArray.size(); i++){
        if(&f == &this->mesh.faceArray[i]){
            index = i;
            break;
        }
    }
    return index;
}

std::list<Face *> ModelAsset::patchToBeDeleted(std::list<int> indexList){
    std::list<Face *> facePatch;
    for(std::list<int>::iterator i = indexList.begin(); i != indexList.end(); i++){
        Face *f = &this->mesh.faceArray[*i];
        facePatch.push_back(f);
    }
    return facePatch;
}

std::list<Face *> ModelAsset::patchToBeDeleted(){
    
    Face *f = &this->mesh.faceArray[0];
    std::list<Face *> facePatch;
    facePatch = f->getFacePatch();
    
//    float qualityThreshold = 10.0f;
//    std::list<Face *> facePatch;
//    Face *faceWithBadQuality = NULL;
//    
//    //Finding one bad quality face...
//    for(int i=0; i<this->mesh.faceArray.size(); i++){
//        Face *f = &this->mesh.faceArray[i];
//        if(f->quality() > qualityThreshold && false == f->processed){
//            //faceWithBadQuality = f;
//            //break;
//            facePatch.push_back(f);
//        }
//    }

//    //finding the initial patch...
//    if(faceWithBadQuality != NULL){
//        std::list<Face *> patch = faceWithBadQuality->getFacePatch();
//        std::list<Face *>::iterator it;
//        
//        for(it = patch.begin(); it != patch.end(); it++){
//            if((**it).quality() > qualityThreshold){
//                facePatch.push_back(*it);
//            }
//        }
//    }
//    
//    //copying the initial patch...
//    std::list<Face *> facePatchMirror = facePatch;
//
//    
//    while(facePatchMirror.size() > 0){
//        Face *f = facePatchMirror.front();
//        std::list<Face *> patch = f->getFacePatch();
//        std::list<Face *>::iterator it;
//        
//        for(it = patch.begin(); it != patch.end(); it++){
//            if((**it).quality() > qualityThreshold){
//                std::list<Face *>::iterator itCheck;
//                bool in = false;
//                for(itCheck = facePatch.begin(); itCheck != facePatch.end(); itCheck++){
//                    if((**it).id == (**itCheck).id){
//                        in = true;
//                        break;
//                    }
//                }
//                if(false == in){
//                    facePatch.push_back(*it);
//                    facePatchMirror.push_back(*it);
//                }
//            }
//        }
//        facePatchMirror.pop_front();
//    }
    return facePatch;
}


Data::Mesh ModelAsset::openHole(){
    
    std::list<Face *> facePatch = this->patchToBeDeleted();
    std::list<Face *>::iterator it = facePatch.begin();
    
    for(it = facePatch.begin(); it != facePatch.end(); it++){
        (**it).toBeErased = true;
        
        int ind = findIndex(**it);

        std::list<Edge>::iterator itEdge = this->mesh.edgeList.begin();
        for(int j=0; j<3*ind; j++){ itEdge++; }
        (*itEdge++).toBeErased = true;
        (*itEdge++).toBeErased = true;
        (*itEdge  ).toBeErased = true;
        
        std::list<Vector>::iterator itNormal = this->mesh.normalList.begin();
        for(int j=0; j<2*ind; j++){ itNormal++; }
        (*itNormal++).toBeErased = true;
        (*itNormal  ).toBeErased = true;
        
        
        for(int i=0; i<3; i++){
            Vector *v = (**it).getVector(i);
            
            std::list<Face *> faceList = v->faceList;
            std::list<Face *>::iterator itFace;
            
            int inPatchCount = 0;
            for(itFace = faceList.begin(); itFace != faceList.end(); itFace++){
                
                bool isInPatch = false;
                
                std::list<Face *>::iterator itPatch;
                for(itPatch = facePatch.begin(); itPatch != facePatch.end(); itPatch++){

                    if(&(**itFace) == &(**itPatch)){
                        isInPatch = true;
                        break;
                    }
                }
                if(isInPatch == true){
                    inPatchCount++;
                }
            }
            
            if(inPatchCount == v->faceList.size()){
                v->toBeErased = true;
            }
        }
    }
    

    Data::Mesh holeMesh;
    holeMesh.faceArray = this->removeToBeErasedFaces();
    //this->removeToBeErasedEdges();
    //this->removeToBeErasedNormals();
    

    holeMesh.makeVertexArrayFromFaceArray();
    
    return holeMesh;
}

Data::Mesh ModelAsset::openHole(std::list<int> indexList){
    
    std::list<Face *> facePatch = this->patchToBeDeleted(indexList);
    std::list<Face *>::iterator it = facePatch.begin();
    
    for(it = facePatch.begin(); it != facePatch.end(); it++){
        (**it).toBeErased = true;
        
        int ind = findIndex(**it);
        
        std::list<Edge>::iterator itEdge = this->mesh.edgeList.begin();
        for(int j=0; j<3*ind; j++){ itEdge++; }
        (*itEdge++).toBeErased = true;
        (*itEdge++).toBeErased = true;
        (*itEdge  ).toBeErased = true;
        
        std::list<Vector>::iterator itNormal = this->mesh.normalList.begin();
        for(int j=0; j<2*ind; j++){ itNormal++; }
        (*itNormal++).toBeErased = true;
        (*itNormal  ).toBeErased = true;
        
        
        for(int i=0; i<3; i++){
            Vector *v = (**it).getVector(i);
            
            std::list<Face *> faceList = v->faceList;
            std::list<Face *>::iterator itFace;
            
            int inPatchCount = 0;
            for(itFace = faceList.begin(); itFace != faceList.end(); itFace++){
                
                bool isInPatch = false;
                
                std::list<Face *>::iterator itPatch;
                for(itPatch = facePatch.begin(); itPatch != facePatch.end(); itPatch++){
                    
                    if(&(**itFace) == &(**itPatch)){
                        isInPatch = true;
                        break;
                    }
                }
                if(isInPatch == true){
                    inPatchCount++;
                }
            }
            
            if(inPatchCount == v->faceList.size()){
                v->toBeErased = true;
            }
        }
    }
    
    Data::Mesh holeMesh;
    holeMesh.faceArray = this->removeToBeErasedFaces();
    //this->removeToBeErasedEdges();
    //this->removeToBeErasedNormals();
    
    //this->removeToBeErasedVertex();
//    for(int i=0; i<this->mesh.faceArray.size(); i++){
//        Face *f = &this->mesh.faceArray[i];
//        Vector *v0 = f->v0;
//        Vector *v1 = f->v1;
//        Vector *v2 = f->v2;
//        
//        int index = -1;
//        for(int j=0; j<this->mesh.vertexArray.size(); j++){
//            Vector *v = &this->mesh.vertexArray[j];
//            if(v->x == v0->x && v->y == v0->y && v->z == v0->z){
//                index = j;
//            }
//        }
//        std::cout << this->mesh.vertexArray[index].toBeErased << std::endl;
//        f->v0 = &this->mesh.vertexArray[index];
//        
//        index = -1;
//        for(int j=0; j<this->mesh.vertexArray.size(); j++){
//            Vector *v = &this->mesh.vertexArray[j];
//            if(v->x == v1->x && v->y == v1->y && v->z == v1->z){
//                index = j;
//            }
//        }
//        std::cout << this->mesh.vertexArray[index].toBeErased << std::endl;
//        f->v1 = &this->mesh.vertexArray[index];
//        
//        index = -1;
//        for(int j=0; j<this->mesh.vertexArray.size(); j++){
//            Vector *v = &this->mesh.vertexArray[j];
//            if(v->x == v2->x && v->y == v2->y && v->z == v2->z){
//                index = j;
//            }
//        }
//        std::cout << this->mesh.vertexArray[index].toBeErased << std::endl;
//        f->v2 = &this->mesh.vertexArray[index];
//    }
    
    holeMesh.makeVertexArrayFromFaceArray();
    return holeMesh;
}

std::vector<Face> ModelAsset::removeToBeErasedFaces(){
    std::vector<Face> faces;
    for(int i=0; i<this->mesh.faceArray.size(); i++){
        if(true == this->mesh.faceArray[i].toBeErased){
            faces.push_back(this->mesh.faceArray[i]);
            this->mesh.faceArray.erase(this->mesh.faceArray.begin() + i);
            i--;
        }
    }
    return faces;
}

void ModelAsset::removeToBeErasedVertex(){
    for(int i=0; i<this->mesh.vertexArray.size(); i++){
        if(true == this->mesh.vertexArray[i].toBeErased){
            this->mesh.vertexArray.erase(this->mesh.vertexArray.begin() + i);
            i--;
        }
    }
}

void ModelAsset::removeToBeErasedEdges(){
    std::list<Edge>::iterator itEdge = this->mesh.edgeList.begin();
    for(itEdge = this->mesh.edgeList.begin(); itEdge != this->mesh.edgeList.end(); itEdge++){
        if( (*itEdge).toBeErased == true ){
            this->mesh.edgeList.erase(itEdge);
        }
    }
}

void ModelAsset::removeToBeErasedNormals(){
    std::list<Vector>::iterator itNormal = this->mesh.normalList.begin();
    for(itNormal = this->mesh.normalList.begin(); itNormal != this->mesh.normalList.end(); itNormal++){
        if( (*itNormal).toBeErased == true ){
            this->mesh.normalList.erase(itNormal);
        }
    }
}

void ModelAsset::putInBuffer(GLuint &vbo, GLuint &vao, GLfloat data[], unsigned long size){
    
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER,vbo);
    glBufferData(GL_ARRAY_BUFFER, size, data, GL_STATIC_DRAW);
    glEnableVertexAttribArray(this->shaders->attrib("vert"));
    glVertexAttribPointer(this->shaders->attrib("vert"), 3, GL_FLOAT, GL_FALSE, 8*sizeof(GLfloat), NULL);
    glEnableVertexAttribArray(this->shaders->attrib("vertTexCoord"));
    glVertexAttribPointer(this->shaders->attrib("vertTexCoord"), 2, GL_FLOAT, GL_TRUE, 8*sizeof(GLfloat), (const GLvoid*)(3 * sizeof(GLfloat)));

    glEnableVertexAttribArray(this->shaders->attrib("vertNormal"));
    glVertexAttribPointer(this->shaders->attrib("vertNormal"), 3, GL_FLOAT, GL_TRUE,  8*sizeof(GLfloat), (const GLvoid*)(5 * sizeof(GLfloat)));

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

int ModelAsset::joinMeshesData(Data::Mesh *newMesh){

    int newVertexAmount = 0;
    
    for(int i=0; i<newMesh->faceArray.size(); i++){
        Face *f = &newMesh->faceArray[i];
        
        for(int j=0; j<3; j++){

            Vector *v = f->getVector(j);
            bool equals = false;
            
            for(int k=0; k<this->mesh.vertexArray.size(); k++){
                Vector *vMesh = &this->mesh.vertexArray[k];
                
                if(v->x == vMesh->x && v->y == vMesh->y && v->z == vMesh->z){
                    equals = true;
                    if(j == 0) f->v0 = vMesh;
                    if(j == 1) f->v1 = vMesh;
                    if(j == 2) f->v2 = vMesh;
                }
            }
            
            //a new vertex has to be pushed into the main mesh vertexArray
            if(false == equals){
                
                int ind = 2;
                Vector last = this->mesh.vertexArray[this->mesh.vertexArray.size()-1];
                while(last.toBeErased == true){
                    last = this->mesh.vertexArray[this->mesh.vertexArray.size()-ind];
                    ind++;
                }
                int id = last.id + 1;
                v->id = id;
                this->mesh.vertexArray.push_back(*v);
                newVertexAmount++;
                
                if(j == 0) f->v0 = &this->mesh.vertexArray[this->mesh.vertexArray.size()-1];
                if(j == 1) f->v1 = &this->mesh.vertexArray[this->mesh.vertexArray.size()-1];
                if(j == 2) f->v2 = &this->mesh.vertexArray[this->mesh.vertexArray.size()-1];
            }
        }
        
        this->mesh.faceArray.push_back(*f);
    }
    
    return newVertexAmount;
    
    
//    this->mesh.faceArray.insert(this->mesh.faceArray.end(),
//                                newMesh->faceArray.begin(),
//                                newMesh->faceArray.end());
    
//    this->mesh.edgeList.insert(this->mesh.edgeList.end(),
//                               newMesh->edgeList.begin(),
//                               newMesh->edgeList.end());
//    
//    this->mesh.normalList.insert(this->mesh.normalList.end(),
//                                 newMesh->normalList.begin(),
//                                 newMesh->normalList.end());
}

void ModelAsset::renderMesh(Data::Mesh *mesh){
    
    int numFaces = static_cast<int>(mesh->faceArray.size());
    GLfloat vertexData[numFaces*3*8];
    
    this->drawCount = numFaces*3;
    
    //for the vertex...
    int k = 0;
    for(int i=0; i<numFaces; i++){
        
        Face *f = &mesh->faceArray[i];
        
        for(int j=0; j<3; j++){
            
            Vector *vec;
            if(j==0) vec = f->v0;
            else if(j==1) vec = f->v1;
            else if(j==2) vec = f->v2;
            
            vertexData[k++] = vec->x;
            vertexData[k++] = vec->y;
            vertexData[k++] = vec->z;
            
            //if(!mesh->getHasTextureCoords()){
                vertexData[k++] = 0.0;
                vertexData[k++] = 0.0;
            //}
            
            vertexData[k++] = f->normal.x;
            vertexData[k++] = f->normal.y;
            vertexData[k++] = f->normal.z;
            
        }
    }
    
//    int     numNormals = static_cast<int>(mesh->normalList.size());
//    GLfloat normalData[numNormals*8];
//
//
//    //for the normals...
//    std::list<Vector>::const_iterator it;
//    it = mesh->normalList.begin();
//    int n = 0;
//    
//    for(int i=0; i<numNormals; i++){
//        
//        normalData[n++] = it->x;
//        normalData[n++] = it->y;
//        normalData[n++] = it->z;
//        
//        normalData[n++] = 0.0;
//        normalData[n++] = 0.0;
//        
//        normalData[n++] = it->x;
//        normalData[n++] = it->y;
//        normalData[n++] = it->z;
//        
//        it++;
//    }
    
    //for the edges...
//    std::list<Edge>::const_iterator itEdge;
//
//    int     numEdges   = static_cast<int>(mesh->edgeList.size());
//    GLfloat edgeData[numEdges*8*2];
//
//    this->drawEdgeCount = numEdges*2;
//
//    int e = 0;
//    for(itEdge = mesh->edgeList.begin(); itEdge != mesh->edgeList.end(); itEdge++){
//        
//        edgeData[e++] = itEdge->v0->x;
//        edgeData[e++] = itEdge->v0->y;
//        edgeData[e++] = itEdge->v0->z;
//        
//        edgeData[e++] = 0.0;
//        edgeData[e++] = 0.0;
//        
//        edgeData[e++] = itEdge->v0->x;
//        edgeData[e++] = itEdge->v0->y;
//        edgeData[e++] = itEdge->v0->z;
//        
//        edgeData[e++] = itEdge->v1->x;
//        edgeData[e++] = itEdge->v1->y;
//        edgeData[e++] = itEdge->v1->z;
//        
//        edgeData[e++] = 0.0;
//        edgeData[e++] = 0.0;
//        
//        edgeData[e++] = itEdge->v1->x;
//        edgeData[e++] = itEdge->v1->y;
//        edgeData[e++] = itEdge->v1->z;
//    }
    
    //VERTEX
    this->putInBuffer(this->vbo, this->vao, vertexData, sizeof(vertexData));
    
    //NORMALS
    //this->putInBuffer(this->normalVbo, this->normalVao, normalData, sizeof(normalData));
    
    //EDGES
    //this->putInBuffer(this->edgeVbo, this->edgeVao, edgeData, sizeof(edgeData));
}

void ModelAsset::renderSelectedFaces(std::list<int> indexList){
    
    std::list<Vector> selectedNormalList;
    for(std::list<int>::iterator i = indexList.begin(); i != indexList.end(); i++){
        Face *f = &this->originalMesh.faceArray[*i];
        selectedNormalList.push_back(f->center());
        Vector n = f->normal.normalize();
        float tam = (*f->v1 - *f->v0).length();
        n.x *= tam;
        n.y *= tam;
        n.z *= tam;
        selectedNormalList.push_back(f->center() + n);
    }
    
    int numNormals = static_cast<int>(indexList.size()*2);
    GLfloat normalData[numNormals*8];

    //for the normals...
    std::list<Vector>::const_iterator it;
    it = selectedNormalList.begin();
    int n = 0;

    for(int i=0; i<numNormals; i++){

        normalData[n++] = it->x;
        normalData[n++] = it->y;
        normalData[n++] = it->z;

        normalData[n++] = 0.0;
        normalData[n++] = 0.0;

        normalData[n++] = it->x;
        normalData[n++] = it->y;
        normalData[n++] = it->z;
        
        it++;
    }
    //NORMALS
    this->putInBuffer(this->normalVbo, this->normalVao, normalData, sizeof(normalData));
    
        
        
    
    std::list<Edge> selectedEdges;
    for(std::list<int>::iterator i = indexList.begin(); i != indexList.end(); i++){
        Face *f = &this->originalMesh.faceArray[*i];
        selectedEdges.push_back(*f->e0);
        selectedEdges.push_back(*f->e1);
        selectedEdges.push_back(*f->e2);
    }
    
    int numEdges = static_cast<int>(indexList.size()*3);
    GLfloat edgeData[numEdges*8*2];
    
    this->drawEdgeCount = numEdges*2;
    
    int e = 0;
    for(std::list<Edge>::iterator itEdge = selectedEdges.begin(); itEdge != selectedEdges.end(); itEdge++){

        edgeData[e++] = itEdge->v0->x;
        edgeData[e++] = itEdge->v0->y;
        edgeData[e++] = itEdge->v0->z;

        edgeData[e++] = 0.0;
        edgeData[e++] = 0.0;

        edgeData[e++] = itEdge->v0->x;
        edgeData[e++] = itEdge->v0->y;
        edgeData[e++] = itEdge->v0->z;

        edgeData[e++] = itEdge->v1->x;
        edgeData[e++] = itEdge->v1->y;
        edgeData[e++] = itEdge->v1->z;

        edgeData[e++] = 0.0;
        edgeData[e++] = 0.0;

        edgeData[e++] = itEdge->v1->x;
        edgeData[e++] = itEdge->v1->y;
        edgeData[e++] = itEdge->v1->z;
    }
    
    
    //EDGES
    //std::cout << sizeof(edgeData ) << std::endl;
    this->putInBuffer(this->edgeVbo, this->edgeVao, edgeData, sizeof(edgeData));

}

void ModelAsset::renderSelectedVertex(std::vector<Vector> initial, std::vector<Vector> end){
    
    std::list<Vector> selectedNormalList;
    
    for(int i=0; i<initial.size(); i++){
        selectedNormalList.push_back(initial[i]);
        selectedNormalList.push_back(end[i]);
    }
    
    int numNormals = static_cast<int>(initial.size()*2);
    GLfloat normalData[numNormals*8];
    
    //for the normals...
    std::list<Vector>::const_iterator it;
    it = selectedNormalList.begin();
    int n = 0;
    
    for(int i=0; i<numNormals; i++){
        
        normalData[n++] = it->x;
        normalData[n++] = it->y;
        normalData[n++] = it->z;
        
        if(i%2 == 0){
            normalData[n++] = 0.0;
            normalData[n++] = 0.0;
        }else{
            normalData[n++] = 1.0;
            normalData[n++] = 1.0;
        }
        
        normalData[n++] = it->x;
        normalData[n++] = it->y;
        normalData[n++] = it->z;
        
        it++;
    }
    //NORMALS
    this->putInBuffer(this->normalVbo, this->normalVao, normalData, sizeof(normalData));
}

tdogl::Program *ModelAsset::LoadShaders(const char *vertFileName, const char *fragFileName){
    std::vector<tdogl::Shader> shadersVec;
    shadersVec.push_back(tdogl::Shader::shaderFromFile(vertFileName, GL_VERTEX_SHADER));
    shadersVec.push_back(tdogl::Shader::shaderFromFile(fragFileName, GL_FRAGMENT_SHADER));
    return new tdogl::Program(shadersVec);
}

void ModelAsset::LoadTexture(const std::string &texturePath){
    tdogl::Bitmap bitmap = tdogl::Bitmap::bitmapFromFile(texturePath);
    bitmap.flipVertically();
    texture = new tdogl::Texture(bitmap);
    textures.push_back(new tdogl::Texture(bitmap));
}

void ModelAsset::LoadCameraMatrix(glm::mat4 &camera){
    this->shaders->use();
    this->camera = &camera;
    this->shaders->setUniform("camera", camera);
    this->shaders->stopUsing();
}

int ModelAsset::getNewVertexAmount(){
    return this->newVertexAmount;
}

Data::Mesh *ModelAsset::getMesh(){
    return &this->mesh;
}
Data::Mesh *ModelAsset::getMshSurfMesh(){
    return &this->mshSurf;
}

Data::Mesh *ModelAsset::getOriginalMesh(){
    return &this->originalMesh;
}

tdogl::Program *ModelAsset::getShaders(){
    return this->shaders;
}

tdogl::Texture *ModelAsset::getTexture(){
    return this->texture;
}
tdogl::Texture *ModelAsset::getEdgeTexture(){
    return this->edgeTexture;
}

glm::mat4 *ModelAsset::getCameraMatrix(){
    return this->camera;
}

glm::mat4 *ModelAsset::getTransformationMatrix(){
    return this->transformation;
}

bool ModelAsset::getLightEnable(){
    return lightEnable;
}

GLuint ModelAsset::getVbo(){
    return this->vbo;
}
GLuint ModelAsset::getVao(){
    return this->vao;
}

GLuint ModelAsset::getNormalVbo(){
    return this->normalVbo;
}

GLuint ModelAsset::getNormalVao(){
    return this->normalVao;
}

GLuint ModelAsset::getEdgeVbo(){
    return this->edgeVbo;
}

GLuint ModelAsset::getEdgeVao(){
    return this->edgeVao;
}

GLuint ModelAsset::getDrawType(){
    return this->drawType;
}
GLuint ModelAsset::getDrawStart(){
    return this->drawStart;
}
GLuint ModelAsset::getDrawCount(){
    return this->drawCount;
}

GLuint ModelAsset::getDrawEdgeCount(){
    return this->drawEdgeCount;
}

GLfloat ModelAsset::getShininess(){
    return this->shininess;
}
glm::vec3 ModelAsset::getSpecularColor(){
    return this->specularColor;
}

void ModelAsset::makeEdgeConnections(Face *f){
    
    f->e0->leftFace = f;
    f->e1->leftFace = f;
    f->e2->leftFace = f;
    
    std::list<Face *> patch = f->getFacePatch();
    
    std::list<Face *>::iterator i;
    
    int howMany = 0;
    
    for(i = patch.begin(); i != patch.end(); i++){
        int cont = 0;
        bool b0 = false;
        bool b1 = false;
        bool b2 = false;
        
        if( (**i).v0 == f->v0 || (**i).v1 == f->v0 || (**i).v2 == f->v0){
            cont++;
            b0 = true;
        }
        
        if( (**i).v0 == f->v1 || (**i).v1 == f->v1 || (**i).v2 == f->v1){
            cont++;
            b1 = true;
        }
        
        if( (**i).v0 == f->v2 || (**i).v1 == f->v2 || (**i).v2 == f->v2){
            cont++;
            b2 = true;
        }
        
        //It found a real Face neighbor
        if(cont == 2) {
            //Edge0
            if(b0 == true && b1 == true){
                f->e0->rightFace = &(**i);
            }
            //Edge1
            else if(b1 == true && b2 == true){
                f->e1->rightFace = &(**i);
            }
            //Edge2
            else if(b2 == true && b0 == true){
                f->e2->rightFace = &(**i);
            }
            
        }
    }
    
    
}


double ModelAsset::phi(std::list<Face *> P){
    
    double max = -1;
    std::list<Face *>::iterator i;
    
    for(i = P.begin(); i != P.end(); i++){
        
        std::list<Face *>::iterator j = i;
        j++;
        std::list<Face *>::iterator k;
        
        if(j != P.end()){
            for(k = j; k != P.end(); k++){
                //glm::vec3 dif = (*k)->normal - (*i)->normal;
                Vector dif = (*k)->normal - (*i)->normal;
                //if(glm::length(dif) > max) max = glm::length(dif);
                if( dif.length() > max) max = dif.length();
            }
        }

    }
    
    return max;
}

double ModelAsset::R(std::list<Face *> P){
    
    double max = -1;
    double sum = 0;
    double epsilon = 0.01;
    
    std::list<Face *>::iterator i;
    for(i = P.begin(); i != P.end(); i++){
        
        std::list<Face *> realNeighborsInsideP = edgeNeighbor(P, **i);
        
        std::list<Face *>::iterator j;
        for(j = realNeighborsInsideP.begin(); j != realNeighborsInsideP.end(); j++){
            
            //glm::vec3 dif = (*i)->normal - (*j)->normal;
            Vector dif = (*i)->normal - (*j)->normal;
            //sum += glm::length(dif);
            sum += dif.length();
            
            //if(glm::length(dif) > max) max = glm::length(dif);
            if(dif.length() > max) max = dif.length();
        }
    }
    
    //necessary because i do every comparison 2 times
    sum /= 2;
    
    return max/(epsilon + sum);
}

std::list<Face *> ModelAsset::edgeNeighbor(std::list<Face *> P,Face &f){
    
    std::list<Face *> realNeighbors;
    std::list<Face *> realNeighborsInsidePatch;
    //std::list<Face *> patch = getFacePatch(f);
    std::list<Face *> patch = f.getFacePatch();
    
    std::list<Face *>::iterator i;
    
    for(i = patch.begin(); i != patch.end(); i++){
        int cont = 0;
        
        if( (**i).v0 == f.v0 || (**i).v1 == f.v0 || (**i).v2 == f.v0){
            cont++;
        }
        
        if( (**i).v0 == f.v1 || (**i).v1 == f.v1 || (**i).v2 == f.v1){
            cont++;
        }
        
        if( (**i).v0 == f.v2 || (**i).v1 == f.v2 || (**i).v2 == f.v2){
            cont++;
        }
        
        if(cont == 2) realNeighbors.push_back(*i);

    }
    
    std::list<Face *>::iterator j;
    for(j = P.begin(); j != P.end(); j++){
        
        std::list<Face *>::iterator k;
        for(k = realNeighbors.begin(); k != realNeighbors.end(); k++){
            if(&(**j) == &(**k)){
                realNeighborsInsidePatch.push_back(*k);
            }
        }
        
    }
    
    return realNeighborsInsidePatch;
}

double ModelAsset::H(std::list<Face *> P){
    return phi(P) * R(P);
}

Vector ModelAsset::patchAverageNormal(std::list<Face *> P){
    
    std::list<Face *>::iterator i;
    //glm::vec3 sumVec;
    Vector sumVec;
    
    for(i = P.begin(); i != P.end(); i++){
        if(i == P.begin()){
            //sumVec  = (**i).normal * (**i).area();
            sumVec.x = (**i).normal.x * (**i).area();
            sumVec.y = (**i).normal.y * (**i).area();
            sumVec.z = (**i).normal.z * (**i).area();
        }
        else{
            //sumVec = sumVec + (**i).normal * (**i).area();
            sumVec.x = sumVec.x + (**i).normal.x * (**i).area();
            sumVec.y = sumVec.y + (**i).normal.y * (**i).area();
            sumVec.z = sumVec.z + (**i).normal.z * (**i).area();
        }
    }

    Vector g;
    g.x = sumVec.x/sumVec.length();
    g.y = sumVec.y/sumVec.length();
    g.z = sumVec.z/sumVec.length();
    return g;
}

void ModelAsset::computeMeshGuidances(){
    
    for (int i=0; i<this->mesh.faceArray.size(); i++){
        
        Face *f = &this->mesh.faceArray[i];
        
        //construct a patch
        std::list<Face *> patch = this->mesh.faceArray[i].getFacePatch();
        
        //compute the consistency function
        double hPatch = H(patch);
        
        //compute the patch average normal
        Vector g = patchAverageNormal(patch);
        
        std::list<Face *>::iterator j;
        for(j = patch.begin(); j != patch.end(); j++){
            if(hPatch < (**j).consistency){
                (**j).consistency = hPatch;
                (**j).guidance = g;
            }
        }

    }
}


void ModelAsset::updatingVertices(){

    
    for (int i = 0; i < this->mesh.faceArray.size(); i++){
        
        //glm::vec3 fNormal = filteredNormal(this->mesh.faceArray[i]);
        Vector fNormal = filteredNormal(this->mesh.faceArray[i]);
        
        Vector *v;
        for(int k = 0; k < 3; k++){
            if(k == 0) v = this->mesh.faceArray[i].v0;
            if(k == 1) v = this->mesh.faceArray[i].v1;
            if(k == 2) v = this->mesh.faceArray[i].v2;
            
            int Fi = v->faceList.size();
            
            //glm::vec3 n;
            Vector n;
            
            std::list<Face *>::iterator it;
            for(it = v->faceList.begin(); it != v->faceList.end(); it++){
                
                Vector dist;
                dist.x = (**it).center().x - v->x;
                dist.y = (**it).center().y - v->y;
                dist.z = (**it).center().z - v->z;
                
                //float num = glm::dot((**it).guidance, dist);
                float num = (**it).guidance.dot(dist);

                if(it == v->faceList.begin()){
                    //n  = (**it).guidance * num;
                    //n = fNormal * num;
                    n.x = fNormal.x * num;
                    n.y = fNormal.y * num;
                    n.z = fNormal.z * num;
                }
                else{
                    //n += (**it).guidance * num;
                    n.x = n.x + fNormal.x * num;
                    n.y = n.y + fNormal.y * num;
                    n.z = n.z + fNormal.z * num;
                }
            }
            
            //n = n/Fi;
            n.x = n.x/Fi;
            n.y = n.y/Fi;
            n.z = n.z/Fi;

            v->x += n.x;
            v->y += n.y;
            v->z += n.z;
        }
    }
    

    
//    Face *f = &this->mesh.faceArray[0];
//    std::list<Face *> facePatch = f->getFacePatch();
//    
//    
//    std::list<Face *>::iterator it = facePatch.begin();
//
//    for(it = facePatch.begin(); it != facePatch.end(); it++){
//        for(int i=0; i<this->mesh.faceArray.size(); i++){
//
//            if( &(**it) == &this->mesh.faceArray[i]){
//                this->mesh.faceArray[i].v0->y += 0.01;
//                this->mesh.faceArray[i].v1->y += 0.01;
//                this->mesh.faceArray[i].v2->y += 0.01;
//            }
//        
//        }
//    }

    
//    std::list<Edge>::iterator bound;
//    for(bound = mesh.boundary.outline.begin(); bound != mesh.boundary.outline.end(); bound++){
//        bound->v0->y += 0.01;
//        bound->v1->y += 0.01;
//    }
    

    
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    GLfloat *data = (GLfloat *) glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);

    int k = 0;
    for(int i = 0; i < this->mesh.faceArray.size(); i++){
        //glm::vec3 normal = this->mesh.faceArray[i].normal;
        Vector normal = this->mesh.faceArray[i].normal;
        data[k]   = this->mesh.faceArray[i].v0->x;
        data[k+1] = this->mesh.faceArray[i].v0->y;
        data[k+2] = this->mesh.faceArray[i].v0->z;
        data[k+5] = normal.x;
        data[k+6] = normal.y;
        data[k+7] = normal.z;
        
        k += 8;
        data[k]   = this->mesh.faceArray[i].v1->x;
        data[k+1] = this->mesh.faceArray[i].v1->y;
        data[k+2] = this->mesh.faceArray[i].v1->z;
        data[k+5] = normal.x;
        data[k+6] = normal.y;
        data[k+7] = normal.z;

        k += 8;
        data[k]   = this->mesh.faceArray[i].v2->x;
        data[k+1] = this->mesh.faceArray[i].v2->y;
        data[k+2] = this->mesh.faceArray[i].v2->z;
        data[k+5] = normal.x;
        data[k+6] = normal.y;
        data[k+7] = normal.z;

        k += 8;
    }
    
    glUnmapBuffer(GL_ARRAY_BUFFER);
}

void ModelAsset::updateEdges(){
    
    glBindBuffer(GL_ARRAY_BUFFER, this->edgeVbo);
    GLfloat *data = (GLfloat *) glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);
    
    int k = 0;
    
    std::list<Edge>::iterator iterator;
    for(iterator = mesh.edgeList.begin(); iterator != mesh.edgeList.end(); iterator++){
        data[k]   = iterator->v0->x;
        data[k+1] = iterator->v0->y;
        data[k+2] = iterator->v0->z;
        k += 8;
        
        data[k]   = iterator->v1->x;
        data[k+1] = iterator->v1->y;
        data[k+2] = iterator->v1->z;
        k += 8;
    }
    
    glUnmapBuffer(GL_ARRAY_BUFFER);
}

void ModelAsset::updateNormals(){
    
    updatingNormalsFromTriangles();
    
    glBindBuffer(GL_ARRAY_BUFFER, this->normalVbo);
    GLfloat *data = (GLfloat *) glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);
    
    int k = 0;
    std::list<Vector>::iterator it;
    for(it = mesh.normalList.begin(); it != mesh.normalList.end(); it++){
        data[k]   = it->x;
        data[k+1] = it->y;
        data[k+2] = it->z;
        k += 8;
    }
    
    glUnmapBuffer(GL_ARRAY_BUFFER);
}

void ModelAsset::updatingNormalsFromTriangles(){
    
    std::list<Vector>::iterator it;
    it = this->mesh.normalList.begin();
    
    for(int i = 0; i < this->mesh.faceArray.size(); i++){
        //updating the normal from new vertex positions
        Vector v01 = *this->mesh.faceArray[i].v1 - *this->mesh.faceArray[i].v0;
        Vector v02 = *this->mesh.faceArray[i].v2 - *this->mesh.faceArray[i].v0;
        
        Vector v01v02 = v01 % v02;
        
        this->mesh.faceArray[i].normal = v01v02.normalize();
//        glm::vec3 c01,c02;
//        c01.x = v01.x;
//        c01.y = v01.y;
//        c01.z = v01.z;
//        
//        c02.x = v02.x;
//        c02.y = v02.y;
//        c02.z = v02.z;
        
        //this->mesh.faceArray[i].normal = glm::normalize(glm::cross(c01,c02));
        
        //glm::vec3 center = this->mesh.faceArray[i].center();
        Vector center = this->mesh.faceArray[i].center();
        //glm::vec3 n = glm::normalize(this->mesh.faceArray[i].normal);
        Vector n = this->mesh.faceArray[i].normal;

        n.x /= NORMAL_FACTOR;
        n.y /= NORMAL_FACTOR;
        n.z /= NORMAL_FACTOR;
        //glm::vec3 tip = this->mesh.faceArray[i].center() + n;
        Vector tip = this->mesh.faceArray[i].center() + n;
        it->x = center.x;
        it->y = center.y;
        it->z = center.z;
        
        it++;
        
        it->x = tip.x;
        it->y = tip.y;
        it->z = tip.z;
        
        it++;

    }
}

Vector ModelAsset::filteredNormal(Face &f){
    //glm::vec3 filteredNormal;
    Vector filteredNormal;
    
    std::list<Face *> patch = f.getFacePatch();
    std::list<Face *>::iterator it;
    
    for(it = patch.begin(); it != patch.end(); it++){
        //glm::vec3 nTemp;
        Vector nTemp;
        float num = (**it).area() * spatialKernel(f.center() , (**it).center()) * rangeKernel(f.guidance, (**it).guidance);
        //nTemp = (**it).normal * num;
        nTemp.x = (**it).normal.x * num;
        nTemp.y = (**it).normal.y * num;
        nTemp.z = (**it).normal.z * num;
        
        if(it == patch.begin())
            filteredNormal = nTemp;
        else
            filteredNormal = filteredNormal + nTemp;
    }
    
    return filteredNormal.normalize();
}

//sigmas between 0.2 and 0.6...
double ModelAsset::spatialKernel(Vector p, Vector q){
    double res;
    double sigmaS = 0.2;
    //res = (glm::length(p-q))/(2 * sigmaS * sigmaS);
    Vector r = p - q;
    res = (r.length())/(2*sigmaS * sigmaS);
    return exp( -res );
}

double ModelAsset::rangeKernel(Vector Ip, Vector Iq){
    double res;
    double sigmaR = 0.2;
    //res = (glm::length(Ip-Iq))/(2 * sigmaR * sigmaR);
    Vector r = Ip - Iq;
    res = (r.length())/(2*sigmaR * sigmaR);

    return exp( -res );
}


















