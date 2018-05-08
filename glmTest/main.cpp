//
//  main.cpp
//  glmTest
//
//  Created by Lucas Andrade on 7/14/15.
//  Copyright (c) 2015 Lucas Andrade. All rights reserved.
//


// third-party libraries
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "glm.hpp"
#include "glm/glm/gtc/matrix_transform.hpp"

// standard C++ libraries
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <time.h>
#include <list>
#include <unistd.h>
#include <sstream>

#include "MainDrive.h"
#include "ModelAsset.h"
#include "Scene.h"

//util classes
#include "MousePicker.h"
#include "ErrorMeasurements.h"

//tdogl classes
#include "Program.h"
#include "Texture.h"
#include "Bitmap.h"
#include "Camera.h"
#include "Light.h"

//some defines
#define GL_POINTS           0x0000
#define GL_LINES            0x0001
#define GL_LINE_LOOP        0x0002
#define GL_LINE_STRIP       0x0003
#define GL_TRIANGLES        0x0004
#define GL_TRIANGLE_STRIP   0x0005
#define GL_TRIANGLE_FAN     0x0006

#define GL_ENABLE_LIGHTINING true


using namespace MeshSurf;

int A = 0;
int B = 1;

//constants
const glm::vec2 SCREEN_SIZE(1024, 768);
const glm::vec2 SHADOWMAPPING_SIZE(1024, 768);
const std::string homePath = "/Users/lucasandrade/Documents/XCode/glmTest/glmTest/";
const std::string finalObjsPath = "/Users/lucasandrade/Documents/XCode/results/objs/";

char *dir = getcwd(NULL, 0);


//globals
GLFWwindow* gWindow = NULL;
tdogl::Camera gCamera;
std::vector<tdogl::Light> gLights;
std::list<asset::ModelAsset *> *gModelInstances = new std::list<asset::ModelAsset *>();

float gDegreesRotated = 0.0;


glm::mat4 cubeTranslate = glm::translate(glm::mat4(), glm::vec3(18,0,18));
glm::mat4 cubeMatrix = glm::scale(cubeTranslate, glm::vec3(4,4,4));
glm::mat4 identity = glm::mat4();


GLuint depthRenderBuffer = 0;
unsigned int renderedTexture;
unsigned int depthTexture;

static glm::vec3 lightInvDir = glm::vec3(-3.37661, 57.2941, 80.9007);
static glm::vec3 lookAtPos = glm::vec3(15,0,15);


glm::mat4 biasMatrix(
                     0.5, 0.0, 0.0, 0.0,
                     0.0, 0.5, 0.0, 0.0,
                     0.0, 0.0, 0.5, 0.0,
                     0.5, 0.5, 0.5, 1.0
                     );


//2d obj without texture
asset::ModelAsset *guiPointer;


unsigned int createTexture(int w, int h, bool isDepth){
    unsigned int textureId;
    glGenTextures(1, &textureId);
    glBindTexture(GL_TEXTURE_2D, textureId);
    glTexImage2D(GL_TEXTURE_2D, 0, (!isDepth ? GL_RGBA8 : GL_DEPTH_COMPONENT32), w, h, 0, isDepth ? GL_DEPTH_COMPONENT : GL_RGBA , GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    
    glBindTexture(GL_TEXTURE_2D, 0);
    
    return textureId;
}

void createFrameBuffer(GLuint &frameBuffer){
    glGenFramebuffers(1, &frameBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);
}


Util::MousePicker *picking;
asset::ModelAsset *modelAssetMesh;
bool stopMouse = false;
bool showNormals = true;
std::list<int> indexList;
bool stop = false;

void UpdateMesh(){
    
    if(glfwGetKey(gWindow,'9') && false == stop){
        stop = true;
        std::cout << "mshsurf" << std::endl;
        modelAssetMesh->executeMshSurf();
    }
    
    if(glfwGetKey(gWindow,'0') && false == stop){
        stop = true;
        if(indexList.size() != 0){
            picking->indexListClear();
            
            /*indexList.clear();
            Data::Mesh* m = modelAssetMesh->getMesh();
            for(int i = 0; i < m->faceArray.size()-1; i++){
                indexList.push_back(i);
            }*/
            std::cout << "teste" << std::endl;
            modelAssetMesh->executeMshSurf(indexList);
            modelAssetMesh->showFinalMesh();
        }
    }
    if(glfwGetKey(gWindow,'1')){
        modelAssetMesh->showFinalMesh();
    }
    if(glfwGetKey(gWindow,'2')){
        modelAssetMesh->showMshSurfMesh();
    }
    if(glfwGetKey(gWindow,'3')){
        modelAssetMesh->showMeshWithHoles();
    }
    if(glfwGetKey(gWindow,'4')){
        modelAssetMesh->showOriginalMesh();
    }
    
    if(glfwGetKey(gWindow,'7')){
        glfwSetCursorPos(gWindow, 0,0);
        stopMouse = false;
    }
    
    if(glfwGetKey(gWindow,'6')){
        glfwSetInputMode(gWindow,GLFW_CURSOR,GLFW_CURSOR_NORMAL);
        glfwSetCursorPos(gWindow, SCREEN_SIZE.x/2,SCREEN_SIZE.y/2);
        stopMouse = true;
    }
    
    
    if(glfwGetMouseButton(gWindow, GLFW_MOUSE_BUTTON_1) && true == stopMouse){
        double mouseX, mouseY;
        glfwGetCursorPos(gWindow, &mouseX, &mouseY);
        int index = picking->update(gCamera.position(), gCamera.view(), mouseX, mouseY);
        picking->addToIndexList(index);
        indexList = picking->getIndexList();
    }
    
    if(glfwGetMouseButton(gWindow, GLFW_MOUSE_BUTTON_2) && true == stopMouse){
        double mouseX, mouseY;
        glfwGetCursorPos(gWindow, &mouseX, &mouseY);
        int index = picking->update(gCamera.position(), gCamera.view(), mouseX, mouseY);
        picking->removeFromIndexList(index);
        indexList = picking->getIndexList();
    }
    
    modelAssetMesh->renderSelectedFaces(indexList);
}

void UpdateLight(){
    
    
    
//    if(glfwGetKey(gWindow,'1')){
//        //gLights[1].setIntensities(glm::vec3(1,1,1));
//        gLights[1].setPosition(glm::vec4(gCamera.position(), 1.0));
//    }
    
    if(glfwGetKey(gWindow,'5')){
        lightInvDir = gCamera.position();
        lookAtPos = gCamera.forward() + gCamera.position();
    }
    
    if(glfwGetMouseButton(gWindow, GLFW_MOUSE_BUTTON_2)){
        //gLights[1].setIntensities(glm::vec3(0,0,0));
        // gLights[2].setIntensities(glm::vec3(0,0,0));
    }
    
    if(glfwGetMouseButton(gWindow, GLFW_MOUSE_BUTTON_1)){
        //gLights[1].setIntensities(glm::vec3(0.3,0.3,0.3));
        //gLights[2].setIntensities(glm::vec3(0.3,0.3,0.3));
        //gLights[2].setPosition(glm::vec4(gCamera.position(), 1.0));
        //gLights[2].setConeDirection(gCamera.forward());
        //showNormals = !showNormals;
    }
    
    if(glfwGetMouseButton(gWindow, GLFW_MOUSE_BUTTON_3)){
        gCamera.setPosition(lightInvDir);
        gCamera.lookAt(lookAtPos);
    }
    
}

void UpdateMouse(){
    glfwSetInputMode(gWindow,GLFW_CURSOR,GLFW_CURSOR_DISABLED);
    const float mouseSensitivity = 0.1f;
    double mouseX, mouseY;
    glfwGetCursorPos(gWindow, &mouseX, &mouseY);
    gCamera.offsetOrientation(mouseSensitivity * (float)mouseY, mouseSensitivity * (float)mouseX);
    glfwSetCursorPos(gWindow, 0,0);
}


float teste = 0;
bool filter = false;
int filterInt = 0;
float moveSpeed = 3.0;

void UpdateKeys(float secondsElapsed){
    
    if(glfwGetKey(gWindow, 'Z')){
        gCamera.offsetPosition(secondsElapsed * moveSpeed * -glm::vec3(0,1,0));
    } else if(glfwGetKey(gWindow, 'X')){
        gCamera.offsetPosition(secondsElapsed * moveSpeed * glm::vec3(0,1,0));
    }
    
    if(glfwGetKey(gWindow, 'W')){
        teste += 0.01;
        gCamera.offsetPosition(secondsElapsed * moveSpeed * glm::cross(glm::vec3(0,1,0), gCamera.right()));
    } if(glfwGetKey(gWindow, 'D')){
        gCamera.offsetPosition(secondsElapsed * moveSpeed * gCamera.right());
    }
    if(glfwGetKey(gWindow, 'S')){
        gCamera.offsetPosition(secondsElapsed * moveSpeed * glm::cross(glm::vec3(0,1,0), -gCamera.right()));
    } if(glfwGetKey(gWindow, 'A')){
        gCamera.offsetPosition(secondsElapsed * moveSpeed * -gCamera.right());
    }
    if(glfwGetKey(gWindow, GLFW_KEY_EQUAL)){
        //filter = !filter;
        moveSpeed += 1.0;
    }
    
    if(glfwGetKey(gWindow, GLFW_KEY_MINUS)){
        if(moveSpeed > 1.0)
            moveSpeed -= 1.0;
    }
    
    if(glfwGetKey(gWindow, GLFW_KEY_TAB )){
    }
}



template <typename T>
void SetLightUniform(tdogl::Program* shaders, const char* propertyName, size_t lightIndex, const T& value) {
    std::ostringstream ss;
    ss << "allLights[" << lightIndex << "]." << propertyName;
    std::string uniformName = ss.str();
    
    shaders->setUniform(uniformName.c_str(), value);
}



void Update(float secondsElapsed) {
    
    const GLfloat degreesPerSecond = 20.0f;
    
    //rotate by 1 degree
    gDegreesRotated += secondsElapsed * degreesPerSecond;
    
    while(gDegreesRotated > 360.0f) gDegreesRotated -= 360.0f;
    
    UpdateKeys(secondsElapsed);
    
    if(!stopMouse) UpdateMouse();
    
    UpdateLight();
    
    UpdateMesh();
    
   // gLights[1].setIntensities(glm::vec3(1,1,1));
    gLights[1].setPosition(glm::vec4(gCamera.position(), 1.0));
}

static void RenderInstance(asset::ModelAsset &asset) {

    tdogl::Program* shaders = asset.getShaders();
    
    //asset.update();


     //Compute the MVP matrix from the light's point of view
    //glm::mat4 depthProjectionMatrix = glm::ortho<float>(-10,10,-10,10,0,10);
    glm::mat4 depthProjectionMatrix = gCamera.projection();
    glm::mat4 depthViewMatrix = glm::lookAt(lightInvDir, lookAtPos, glm::vec3(0,1,0));
    glm::mat4 depthModelMatrix = *asset.getTransformationMatrix();
    glm::mat4 depthMVP = depthProjectionMatrix * depthViewMatrix * depthModelMatrix;
    glm::mat4 depthBiasMVP = biasMatrix*depthMVP;

    

    if(filter){
        filter = !filter;
        asset.computeMeshGuidances();
        asset.updatingVertices();
        asset.updateNormals();
        asset.updateEdges();
    }
    
//    if(filterInt < 200){
//        asset.computeMeshGuidances();
//        asset.updatingVertices();
//        asset.updateNormals();
//        asset.updateEdges();
//        filterInt++;
//    }
    
    //bind the shaders
    shaders->use();
    
    if(asset.getTexture() != NULL){
        // bind the texture and set the "tex" uniform in the fragment shader
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, asset.getTexture()->object());
        shaders->setUniform("materialTex", 0); //set to 0 because the texture is bound to GL_TEXTURE0
    }
    
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, depthTexture);
    shaders->setUniform("shadowMap", 1);
    
    //apply the asset transformation
    if(asset.getTransformationMatrix() != NULL){
        shaders->setUniform("model", *asset.getTransformationMatrix());
    }
    
    //apply camera transformations
    if(asset.getCameraMatrix() != NULL){
        
        shaders->setUniform("camera", gCamera.matrix());
        shaders->setUniform("cameraPosition", gCamera.position());
    }
    shaders->setUniform("DepthBiasMVP", depthBiasMVP);
    
    
    //enable lightning on shader
    if(asset.getLightEnable()){
        //heavy calculations in fragment for each pixel, so it is better to compute it here
        glm::mat3 normalMatrix = glm::transpose(glm::inverse(glm::mat3(*asset.getTransformationMatrix())));
        shaders->setUniform("normalMatrix", normalMatrix);
        
        shaders->setUniform("numLights", (int)gLights.size());
        
        for(size_t i = 0; i < gLights.size(); ++i){
            SetLightUniform(shaders, "position", i, gLights[i].getPosiion());
            SetLightUniform(shaders, "intensities", i, gLights[i].getIntensities());
            SetLightUniform(shaders, "attenuation", i, gLights[i].getAttenuation());
            SetLightUniform(shaders, "ambientCoefficient", i, gLights[i].getAmbientCoefficient());
            SetLightUniform(shaders, "coneAngle", i, gLights[i].getConeAngle());
            SetLightUniform(shaders, "coneDirection", i, gLights[i].getConeDirection());
        }
        
        //material related
        shaders->setUniform("materialShininess", asset.getShininess());
        shaders->setUniform("materialSpecularColor", asset.getSpecularColor());
    }
    
    //bind VAO and draw
    glBindVertexArray(asset.getVao());
    glDrawArrays(asset.getDrawType(), asset.getDrawStart(), asset.getDrawCount());
    
    
    if(asset.getEdgeTexture() != NULL){
        // bind the texture and set the "tex" uniform in the fragment shader
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, asset.getEdgeTexture()->object());
        shaders->setUniform("materialTex", 1); //set to 0 because the texture is bound to GL_TEXTURE0
    }
    
    
   if(showNormals && indexList.size() != 0){
        //showNormals = !showNormals;
        glBindVertexArray(asset.getNormalVao());
        glDrawArrays(GL_LINES, asset.getDrawStart(), asset.getDrawCount());
    }
    

   glBindVertexArray(asset.getEdgeVao());
    glDrawArrays(GL_LINES, asset.getDrawStart(), asset.getDrawEdgeCount());

    
    //unbind everything
    glBindVertexArray(0);
    shaders->stopUsing();
    
}



static void RenderShadowmap(asset::ModelAsset &asset){
    
    std::vector<tdogl::Shader> shadersVec;
    shadersVec.push_back(tdogl::Shader::shaderFromFile(homePath + "shaders/shadowmapping2.vs", GL_VERTEX_SHADER));
    shadersVec.push_back(tdogl::Shader::shaderFromFile(homePath + "shaders/shadowmapping2.fs", GL_FRAGMENT_SHADER));
    tdogl::Program* shaders = new tdogl::Program(shadersVec);
    
    // Compute the MVP matrix from the light's point of view
    //glm::mat4 depthProjectionMatrix = glm::ortho<float>(-10,10,-10,10,0,10);
    glm::mat4 depthProjectionMatrix = gCamera.projection();
    glm::mat4 depthViewMatrix = glm::lookAt(lightInvDir, lookAtPos, glm::vec3(0,1,0));
    glm::mat4 depthModelMatrix = *asset.getTransformationMatrix();
    glm::mat4 depthMVP = depthProjectionMatrix * depthViewMatrix * depthModelMatrix ;
    
    
    shaders->use();
    
    shaders->setUniform("depthMVP",depthMVP);
    //shaders->setUniform("depthMVP", gCamera.matrix() * depthModelMatrix);


    
    glBindVertexArray(asset.getVao());
    glDrawArrays(asset.getDrawType(), asset.getDrawStart(), asset.getDrawCount());
    
    //unbind everything
    glBindVertexArray(0);    
    shaders->stopUsing();
}

static void RenderGui(asset::ModelAsset *asset ) {
    
    tdogl::Program* shaders = asset->getShaders();

    //bind the shaders
    shaders->use();
    
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, depthTexture);
    shaders->setUniform("texture2D", 1);
    
    //bind VAO and draw
    glBindVertexArray(asset->getVao());
    glDrawArrays(asset->getDrawType(), asset->getDrawStart(), asset->getDrawCount());
    
    //unbind everything
    glBindVertexArray(0);
    shaders->stopUsing();
    
}


void Render(){
    
    gLights[0].setPosition(glm::vec4(lightInvDir.x,lightInvDir.y,lightInvDir.z,1));
    
    // render all the instances
    std::list<asset::ModelAsset *>::iterator it;
    
    
    //render on the frame buffer
    glBindFramebuffer(GL_FRAMEBUFFER, depthRenderBuffer);
    glViewport(0,0,SHADOWMAPPING_SIZE.x, SHADOWMAPPING_SIZE.y);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |  GL_STENCIL_BUFFER_BIT);
    

    for(it = gModelInstances->begin(); it != gModelInstances->end(); ++it){
        RenderShadowmap(**it);
    }
    
    //render on screen
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0,0,SCREEN_SIZE.x,SCREEN_SIZE.y);
    //glEnable(GL_CULL_FACE);
    //glCullFace(GL_BACK); // Cull back-facing triangles -> draw only front-facing triangles
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    

    for(it = gModelInstances->begin(); it != gModelInstances->end(); ++it){
        RenderInstance(**it);
    }
    
    //RenderGui(guiPointer);

    
    // swap the display buffers (displays what was just drawn)
    glfwSwapBuffers(gWindow);
    glfwPollEvents();
    
}



void OnError(int errorCode, const char* msg) {
    throw std::runtime_error(msg);
}

void InitializeCamera(){
    gCamera.setPosition(glm::vec3(6,2,0));
    gCamera.offsetOrientation(20, 150);
    gCamera.setViewportAspectRatio(SCREEN_SIZE.x / SCREEN_SIZE.y);
    gCamera.setNearAndFarPlanes(0.001, 500);
    gCamera.lookAt(glm::vec3(0,0,0));
}

void InitializeLight(){
    
    //setup lights
    tdogl::Light spotlight;
    spotlight.setPosition(glm::vec4(-4,0,10,1));
    spotlight.setIntensities(glm::vec3(0.1,0.1,0.1)); //strong white light
    spotlight.setAttenuation(0.0001f);
    spotlight.setAmbientCoefficient(0.0f); //no ambient light
    spotlight.setConeAngle(20.0f);
    spotlight.setConeDirection(glm::vec3(0,0,-1));
    
    tdogl::Light directionalLight1;
    directionalLight1.setPosition(glm::vec4(-4,0,10,0)); //w == 0 indications a directional light
    directionalLight1.setIntensities(glm::vec3(0.6,0.6,0.6));
    directionalLight1.setAttenuation(0.0001f);
    directionalLight1.setAmbientCoefficient(0.009f);
    
    tdogl::Light directionalLightShadowMapping;
    directionalLightShadowMapping.setPosition(glm::vec4(lightInvDir, 0)); //w == 0 indications a directional light
    directionalLightShadowMapping.setIntensities(glm::vec3(0.1,0.1,0.1));
    directionalLightShadowMapping.setAttenuation(0.0001f);
    directionalLightShadowMapping.setAmbientCoefficient(0.0009f);

    gLights.push_back(directionalLightShadowMapping);
    gLights.push_back(directionalLight1);
}

void InitializeGLFWandGLEW(){
    
    // initialise GLFW
    glfwSetErrorCallback(OnError);
    if(!glfwInit())
        throw std::runtime_error("glfwInit failed");


    //open a window with GLFW
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    //gWindow = glfwCreateWindow(SCREEN_SIZE.x, SCREEN_SIZE.y, "OpenGL", glfwGetPrimaryMonitor(), NULL);
    gWindow = glfwCreateWindow(SCREEN_SIZE.x, SCREEN_SIZE.y, "OpenGL", NULL, NULL);
    
    // GLFW settings
    glfwSetInputMode(gWindow, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetCursorPos(gWindow, 0, 0);
    
    if(!gWindow)
        throw std::runtime_error("glfwCreateWindow failed. Can your hardware handle OpenGL 3.2?");

    //GLFW settings
    glfwMakeContextCurrent(gWindow);

    //initialise GLEW
    glewExperimental = GL_TRUE; //stops glew crashing on OSX :-/
    if(glewInit() != GLEW_OK)
        throw std::runtime_error("glewInit failed");

    //make sure OpenGL version 3.2 API is available
    if(!GLEW_VERSION_3_2)
        throw std::runtime_error("OpenGL 3.2 API is not available.");
}

void InitializeParameters(){
    
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

}

int main(int argc, const char * argv[]) {
    
    double lastTime = glfwGetTime();
    double thisTime;
    
    InitializeGLFWandGLEW();
    InitializeParameters();
    InitializeCamera();
    InitializeLight();
    
//    std::vector<asset::ModelAsset> vAssets;
//
//    asset::ModelAsset gui(homePath + "objs/PointsFile.txt",
//                          homePath + "shaders/Passthrough.vs",
//                          homePath + "shaders/WobblyTexture.fs",
//                          GL_TRIANGLES);
//    
//    guiPointer = &gui;
//
//    asset::Scene scene(homePath + "objs/secondScene.obj",
//                       homePath + "shaders/3d-lighted-textured-vertex-shader.vs",
//                       homePath + "shaders/3d-lighted-textured-fragment-shader.fs",
//                       GL_TRIANGLES,
//                       GL_ENABLE_LIGHTINING,
//                       gCamera.matrix(),
//                       cubeMatrix);
//
//    vAssets = scene.getAssets();
//    for(int i = 0; i<vAssets.size(); i++){
//        gModelInstances->push_back(vAssets.at(i));
//    }
//
    
    


    
//    asset::ModelAsset obj(//homePath + "objs/test/simplex/simplex.obj",
//                          homePath + "objs/test/block4/filtered_with_hole.obj",
//                          homePath + "shaders/3d-lighted-textured-vertex-shader.vs",
//                          homePath + "shaders/3d-lighted-textured-fragment-shader.fs",
//                          homePath + "textures/yellow.jpg",
//                          homePath + "textures/red.png",
//                          GL_TRIANGLES,
//                          GL_ENABLE_LIGHTINING,
//                          gCamera.matrix(),
//                          identity,
//                          70.0,
//                          glm::vec3(0.5,0.5,0.5));
//    
//    obj.readBaseMesh(homePath + "objs/test/block4/hole_noised.obj");
//    modelAssetMesh = &obj;
//    gModelInstances->push_back(&obj);
    

    
    //std::string filesPath = "/Users/lucasandrade/Documents/XCode/glmTest/glmTest/";

    //finais/mechanic/mechanic_pre_filtered.obj
//    asset::ModelAsset obj(finalObjsPath + "finais/mechanic/mechanic_fast_filtered.obj",
//                          homePath + "shaders/3d-lighted-textured-vertex-shader.vs",
//                          homePath + "shaders/3d-lighted-textured-fragment-shader.fs",
//                          homePath + "textures/yellow.jpg",
//                          homePath + "textures/red2.png",
//                          GL_TRIANGLES,
//                          GL_ENABLE_LIGHTINING,
//                          gCamera.matrix(),
//                          identity,
//                          70.0,
//                          glm::vec3(0.5,0.5,0.5));
//
//    modelAssetMesh = &obj;
//    gModelInstances->push_back(&obj);
    
    
        asset::ModelAsset obj(finalObjsPath + "new_simplex.obj",
                              homePath + "shaders/3d-lighted-textured-vertex-shader.vs",
                              homePath + "shaders/3d-lighted-textured-fragment-shader.fs",
                              homePath + "textures/yellow.jpg",
                              homePath + "textures/red2.png",
                              GL_TRIANGLES,
                              GL_ENABLE_LIGHTINING,
                              gCamera.matrix(),
                              identity,
                              70.0,
                              glm::vec3(0.5,0.5,0.5));
    
        modelAssetMesh = &obj;
        gModelInstances->push_back(&obj);

    

//    asset::ModelAsset groundTruth(finalObjsPath + "finais/mechanic/mechanic_gt.obj",
//                                  homePath + "shaders/3d-lighted-textured-vertex-shader.vs",
//                                  homePath + "shaders/3d-lighted-textured-fragment-shader.fs",
//                                  homePath + "textures/blue.png",
//                                  homePath + "textures/red.png",
//                                  GL_TRIANGLES,
//                                  GL_ENABLE_LIGHTINING,
//                                  gCamera.matrix(),
//                                  identity,
//                                  70.0,
//                                  glm::vec3(0.5,0.5,0.5));
//    gModelInstances->push_back(&groundTruth);
//    
//    
//    Util::ErrorMeasurements error(groundTruth, obj);
//    std::cout << error.volumeError() << std::endl;
    

//    Util::ErrorMeasurements errorM(groundTruth, obj, homePath + "objs/test/block7/aux.txt", homePath + "objs/test/block7/hole.txt");
//    double error = errorM.meanDistanceError();
//    //double error = errorM.meanDistanceToFaces(homePath + "objs/test/block7/ground_truth.obj", true);
//    std::cout << "mean error: " << error << "\n" << std::endl;
//
//    modelAssetMesh->renderSelectedVertex(errorM.initial, errorM.end);
    
    
    
    //Picking 3d...
    asset::ModelAsset *asset = modelAssetMesh;
    Data::Mesh *mesh = asset->getOriginalMesh();
    glm::mat4 projection = gCamera.projection();
    Util::MousePicker mousePicker(projection, static_cast<float>(SCREEN_SIZE.x), static_cast<float>(SCREEN_SIZE.y), mesh);
    picking = &mousePicker;
    //Picking 3d...

    depthTexture = createTexture(SHADOWMAPPING_SIZE.x,SHADOWMAPPING_SIZE.y, true);
    renderedTexture = createTexture(SHADOWMAPPING_SIZE.x,SHADOWMAPPING_SIZE.y, false);
    

    createFrameBuffer(depthRenderBuffer);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthTexture,0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, renderedTexture, 0);

    
    glEnable(GL_LINE_WIDTH);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(10.0f);
    
    glClearColor(1,1,1,0);
    // run while the window is open
    while(!glfwWindowShouldClose(gWindow)){
        // process pending events
        glfwPollEvents();

        // changes to be made before the Render function
        thisTime = glfwGetTime();
        Update((float) (thisTime - lastTime));
        lastTime = thisTime;
        
        // draw one frame
        Render();
        
        
        if(glfwGetKey(gWindow, GLFW_KEY_ESCAPE))
            glfwSetWindowShouldClose(gWindow, GL_TRUE);
    }
    
    // clean up and exit
    glfwTerminate();

    return 0;
}
