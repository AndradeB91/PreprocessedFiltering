//
//  Scene.cpp
//  glmTest
//
//  Created by Lucas Andrade on 1/18/16.
//  Copyright (c) 2016 Lucas Andrade. All rights reserved.
//

#include "Scene.h"

#define SHININESS_DEFAULT 60

using namespace asset;
using namespace std;

Scene::Scene(const std::string &objFilePath,
             const std::string &vertexShader,
             const std::string &fragmentShader,
                const GLenum drawType,
             const bool light,
             glm::mat4 camera,
             glm::mat4 &transformation){
    
    Assimp::Importer importer;
    const char *file_path = strcpy((char*)malloc(objFilePath.length()+1), objFilePath.c_str());
    const aiScene *scene = importer.ReadFile(file_path,0);
    string texturePath = "/Users/lucasandrade/Documents/XCode/glmTest/glmTest/textures/";
    
    for(int i=0; i<scene->mNumMeshes; i++){
        aiMesh *mesh = scene->mMeshes[i];
        const aiMaterial *material = scene->mMaterials[mesh->mMaterialIndex];
        int texIndex = 0;
        aiString path;
        string sTextureName;
        glm::vec3 ambColor, difColor, spColor;
        float shininess;
        
        if(AI_SUCCESS == material->GetTexture(aiTextureType_DIFFUSE, texIndex, &path) ){
            sTextureName = texturePath + path.data;
        }
        
        aiColor4D ambient;
        if(AI_SUCCESS == aiGetMaterialColor(material , AI_MATKEY_COLOR_SPECULAR, &ambient)){
            ambColor.x = ambient.r;
            ambColor.y = ambient.g;
            ambColor.z = ambient.b;
        }
        aiColor4D diffuse;
        if(AI_SUCCESS == aiGetMaterialColor(material , AI_MATKEY_COLOR_SPECULAR, &diffuse)){
            difColor.x = diffuse.r;
            difColor.y = diffuse.g;
            difColor.z = diffuse.b;
        }
        aiColor4D specular;
        if(AI_SUCCESS == aiGetMaterialColor(material , AI_MATKEY_COLOR_SPECULAR, &specular)){
            spColor.x = specular.r;
            spColor.y = specular.g;
            spColor.z = specular.b;
        }
        
        if(0 != material->Get(AI_MATKEY_SHININESS, shininess)){
            shininess /= 4.f;
        }
        else{
            shininess = SHININESS_DEFAULT;
        }
        
        
        ModelAsset m(mesh,vertexShader,fragmentShader,sTextureName,drawType,light,camera,transformation,shininess,spColor,mesh->HasTextureCoords(0));
        assets.push_back(m);
    }
}

std::vector<ModelAsset> Scene::getAssets(){
    return this->assets;
}
