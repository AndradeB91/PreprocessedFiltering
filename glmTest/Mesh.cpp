//
//  Mesh.cpp
//  glmTest
//
//  Created by Lucas Andrade on 30/01/17.
//  Copyright © 2017 Lucas Andrade. All rights reserved.
//

#include "Mesh.h"

using namespace Data;

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


Mesh::Mesh(){
    this->name = "/Users/lucasandrade/Documents/XCode/glmTest/glmTest/objs/default.obj";
}

Mesh::~Mesh(){
    
}

void Mesh::setName(std::string name){
    this->name = name;
}


void Mesh::findOutline(){
    for(int i=0; i<this->faceArray.size(); i++){
        Face *f = &this->faceArray[i];
        
        if(f->e0->rightFace == NULL){
            this->boundary.outline.push_back(*f->e0);
        }
        
        if(f->e1->rightFace == NULL){
            this->boundary.outline.push_back(*f->e1);
        }
        
        if(f->e2->rightFace == NULL){
            this->boundary.outline.push_back(*f->e2);
        }
    }
}

int Mesh::at(aiVector3D vec){
    double EPSILON = 0.00001;
        for(int i=0; i<this->vertexArray.size(); i++){
        Vector *v = &this->vertexArray[i];
        
        if(fabs(vec.x - v->x) <= EPSILON &&
           fabs(vec.y - v->y) <= EPSILON &&
           fabs(vec.z - v->z) <= EPSILON){
            return i;
        }
    }
    return -1;
}

int Mesh::at(Vector vec){
    for(int i=0; i<this->vertexArray.size(); i++){
        if(vec.x == this->vertexArray[i].x &&
           vec.y == this->vertexArray[i].y &&
           vec.z == this->vertexArray[i].z){
            return i;
        }
    }
    return -1;
}

void Mesh::makeEdgeConnections(Face *f){
    
    f->e0->leftFace = f;
    f->e1->leftFace = f;
    f->e2->leftFace = f;
    
    std::list<Face *> patch = f->getFacePatch();
    
    std::list<Face *>::iterator i;
    
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

float Mesh::volume(){
    float volume = 0.0;
    for(int i=0; i<this->faceArray.size(); i++){
        Face *f = &this->faceArray[i];
        volume += f->signedVolume();
    }
    return volume;
}

int Mesh::countReps(Vector v){
    int count = 0;
    for(int i=0; i<this->vertexArray.size(); i++){
        Vector vec = this->vertexArray[i];
        if(v.x == vec.x && v.y == vec.y && v.z == vec.z){
            count++;
        }
    }
    return count;
}

void Mesh::remakeFacePointers(){
    for(int i=0; i<this->faceArray.size(); i++){
        Face *f = &this->faceArray[i];
        for(int j=0; j<3; j++){
            Vector v = *f->getVector(j);
            int index = at(v);
            if     (j == 0) f->v0 = &this->vertexArray[index];
            else if(j == 1) f->v1 = &this->vertexArray[index];
            else if(j == 2) f->v2 = &this->vertexArray[index];
        }
    }
}

void Mesh::makeVertexArrayFromDuplicatedFaceArray(){
    
//    this->vertexArray.clear();
//    for(int i=0; i<this->faceArray.size(); i++){
//        for(int j=0; j<3; j++){
//            Vector *v = this->faceArray[i].getVector(j);
//            if(j == 0)       this->vertexArray.push_back(*this->faceArray[i].v0);
//            else if(j == 1)  this->vertexArray.push_back(*this->faceArray[i].v1);
//            else if(j == 2)  this->vertexArray.push_back(*this->faceArray[i].v2);
//            
//            //this->vertexArray.push_back(*v);
//            std::cout << v->x << " " << v->y << " " << v->z << std::endl;
//        }
//    }
//    this->vertexArray.clear();
//    for(int i=0; i<this->faceArray.size(); i++){
//        for(int j=0; j<3; j++){
//            Vector v = *this->faceArray[i].getVector(j);
//            this->vertexArray.push_back(v);
//            std::cout << v.x << " "<< v.y << " " << v.z << std::endl;
//
//        }
//    }
//
//    std::vector<Vector> vertexArrayTemp = this->vertexArray;
//    //this->vertexArray.clear();
//    for(int i=0; i<vertexArrayTemp.size(); i++){
//        Vector v = this->vertexArray[i];
//        //std::cout <<"reps: " << this->countReps(v) << " " << v.x << " "<< v.y << " " << v.z << std::endl;
////        Vector v = vertexArrayTemp[i];
////        if(!this->isInVertexArray(vertexArrayTemp[i])){
////            this->vertexArray.push_back(v);
////        }
//    }
//    
//    std::cout << " teste: " << this->vertexArray.size() << std::endl;
//    
////    for(int i=0; i<this->vertexArray.size(); i++){
////        int cont = 0;
////        Vector v = this->vertexArray[i];
////        for(int j=0; j<this->vertexArray.size(); j++){
////            Vector vec = this->vertexArray[j];
////            if(v.x == vec.x && v.y == vec.y && v.z == vec.z){
////                cont++;
////                std::cout << v.x << " " << v.y << " " << v.z << std::endl;
////            }
////        }
////        if(cont > 1){
////            this->vertexArray.erase(this->vertexArray.begin()+i);
////            i--;
////        }
////    }
    
}
void Mesh::makeVertexArrayFromFaceArray(){
    this->vertexArray.clear();
    for(int i=0; i<this->faceArray.size(); i++){
        for(int j=0; j<3; j++){
            Vector v = *this->faceArray[i].getVector(j);
            if(!this->isInVertexArray(v)){
                this->vertexArray.push_back(v);
            }
        }
    }
}

bool Mesh::isInVertexArray(Vector v){
    for(int i=0; i<this->vertexArray.size(); i++){
        Vector vec = this->vertexArray[i];
        if(v.x == vec.x && v.y == vec.y && v.z == vec.z){
            return true;
        }
    }
    return false;
}

bool Mesh::constructDataStructureForHoleMesh(){
    
    for(int i=0; i<this->faceArray.size(); i++){
        Face *f = &this->faceArray[i];
        f->e0->leftFace = NULL;
        f->e0->rightFace = NULL;
        f->e1->leftFace = NULL;
        f->e1->rightFace = NULL;
        f->e2->leftFace = NULL;
        f->e2->rightFace = NULL;
    }
    

    
    this->hasNormals = false;
    this->hasTextureCoords = false;
    
    for(int i=0; i<this->faceArray.size(); i++){
        Face *f = &this->faceArray[i];
        
        f->consistency = INT_MAX;
        this->normalList.push_back(f->center());
        this->normalList.push_back(f->center()+f->normal);
    }
    
    for(int i=0; i<this->faceArray.size(); i++){
        Face *f = &this->faceArray[i];
        
        Vector *v0 = f->v0;
        Vector *v1 = f->v1;
        Vector *v2 = f->v2;
        
        v0->faceList.clear();
        v1->faceList.clear();
        v2->faceList.clear();
    }
    //Constructing the backwards pointers (v->fi)
    for(int i=0; i<this->faceArray.size(); i++){
        
        Face *f = &this->faceArray[i];
        
        Vector *v0 = f->v0;
        Vector *v1 = f->v1;
        Vector *v2 = f->v2;
        
        v0->faceList.push_back(f);
        v1->faceList.push_back(f);
        v2->faceList.push_back(f);
    }
    
    //Constructing the edges and the face to edge pointers (f->e = e)
    for(int i=0; i<this->faceArray.size(); i++){
        
        Face *f = &this->faceArray[i];
        
        Edge *e0 = new Edge(f->v0, f->v1);
        Edge *e1 = new Edge(f->v1, f->v2);
        Edge *e2 = new Edge(f->v2, f->v0);
        
        f->e0 = e0;
        f->e1 = e1;
        f->e2 = e2;
        
        //adding the face edges to the mesh...
        this->edgeList.push_back(*e0);
        this->edgeList.push_back(*e1);
        this->edgeList.push_back(*e2);
    }
    
    //Constructing edge to faces pointers (e->leftFace = f1 && e->rightFace = f2)
    for(int i=0; i<this->faceArray.size(); i++){
        
        Face *f = &this->faceArray[i];
        this->makeEdgeConnections(f);
    }
    
    //Finding the border edges of the mesh
    this->findOutline();
    
    return true;
}

bool Mesh::convertFromAssimpToMesh_simpleVersion(aiMesh *assimpMesh){
    int idGen = 0;
    this->hasNormals = false;
    this->hasTextureCoords = false;
    
    for(int i=0; i<assimpMesh->mNumVertices; i++){
        Vector vec;
        vec.x = assimpMesh->mVertices[i].x;
        vec.y = assimpMesh->mVertices[i].y;
        vec.z = assimpMesh->mVertices[i].z;
        vec.id = idGen++;
        this->vertexArray.push_back(vec);
    }
    
    
    int id = 0;
    for(unsigned int i=0; i<assimpMesh->mNumFaces; i++){
        const aiFace& face = assimpMesh->mFaces[i];
        
        Vector faceNormal;
        if(!this->hasNormals){
            aiVector3D v01 = assimpMesh->mVertices[face.mIndices[1]] - assimpMesh->mVertices[face.mIndices[0]];
            aiVector3D v02 = assimpMesh->mVertices[face.mIndices[2]] - assimpMesh->mVertices[face.mIndices[0]];
            Vector c01, c02, c01c02;
            c01.x = v01.x;
            c01.y = v01.y;
            c01.z = v01.z;
            
            c02.x = v02.x;
            c02.y = v02.y;
            c02.z = v02.z;
            
            c01c02 = c01 % c02;
            faceNormal = c01c02.normalize();
        }
        
        Face f;
        int index;
        for(int j=0; j<3; j++){
            
            aiVector3D pos = assimpMesh->mVertices[face.mIndices[j]];
            
            Vector *vec = &this->vertexArray[face.mIndices[j]];
            
            if     (j == 0) f.v0 = vec;
            else if(j == 1) f.v1 = vec;
            else if(j == 2) f.v2 = vec;
            
        }
    
        f.normal = faceNormal;
        f.consistency = INT_MAX;
        f.id = id++;
        this->faceArray.push_back(f);
    }
    
    return true;
}


bool Mesh::convertFromAssimpToMesh(aiMesh *assimpMesh, std::string objFilePath){
    
    int idGen = 0;
    this->hasNormals = assimpMesh->HasNormals();
    this->hasTextureCoords = assimpMesh->HasTextureCoords(0);

    this->hasNormals = false;
    this->hasTextureCoords = false;
    
    
    std::string data, x,y,z;
    std::ifstream file(objFilePath);
    
    while (file >> data) {

        if(data == "v"){
            file >> x;
            file >> y;
            file >> z;
            
            float xValue = std::atof(x.c_str());
            float yValue = std::atof(y.c_str());
            float zValue = std::atof(z.c_str());
            
            Vector vec(xValue, yValue, zValue);
            vec.id = idGen++;
            this->vertexArray.push_back(vec);
            //std::cout << vec.x << " " << vec.y << " " << vec.z << std::endl;

        }
    }
    
    
    
//    //putting all the aiVector3D in vertexArray
//    for(int i=0; i<assimpMesh->mNumVertices; i++){
//        Vector vec;
//        vec.x = assimpMesh->mVertices[i].x;
//        vec.y = assimpMesh->mVertices[i].y;
//        vec.z = assimpMesh->mVertices[i].z;
//        
//        //if(!this->isInVertexArray(vec)){
//            vec.id = idGen++;
//            this->vertexArray.push_back(vec);
//        //}
//    }
    

    int id = 0;
    for(unsigned int i=0; i<assimpMesh->mNumFaces; i++){
        const aiFace& face = assimpMesh->mFaces[i];
        
        Vector faceNormal;
        if(!this->hasNormals){
            aiVector3D v01 = assimpMesh->mVertices[face.mIndices[1]] - assimpMesh->mVertices[face.mIndices[0]];
            aiVector3D v02 = assimpMesh->mVertices[face.mIndices[2]] - assimpMesh->mVertices[face.mIndices[0]];
            Vector c01, c02, c01c02;
            c01.x = v01.x;
            c01.y = v01.y;
            c01.z = v01.z;
            
            c02.x = v02.x;
            c02.y = v02.y;
            c02.z = v02.z;
            
            c01c02 = c01 % c02;
            faceNormal = c01c02.normalize();
        }
        
        Face f;
        int index;
        
        for(int j=0; j<3; j++){
            
            aiVector3D pos = assimpMesh->mVertices[face.mIndices[j]];
            index = at(pos);
            
            Vector *vec = &this->vertexArray[index];
    
            
            if     (j == 0) f.v0 = vec;
            else if(j == 1) f.v1 = vec;
            else if(j == 2) f.v2 = vec;
            
            
            if(!this->hasNormals){
                f.normal = faceNormal;
                
            }else{
                aiVector3D normal = assimpMesh->mNormals[face.mIndices[j]];
                Vector faceNormal;

                
                faceNormal.x = normal.x;
                faceNormal.y = normal.y;
                faceNormal.z = normal.z;
                f.normal = faceNormal;
            }
        }
        
        this->normalList.push_back(f.center());
        Vector n = faceNormal.normalize();
        this->normalList.push_back(f.center() + n);
        
        f.consistency = INT_MAX;
        f.id = id++;
        this->faceArray.push_back(f);
    }
    
    std::cout << "teste" << std::endl;

    //Constructing the backwards pointers (v->fi)
    for(int i=0; i<this->faceArray.size(); i++){
        

        Face *f = &this->faceArray[i];
        
//        std::cout << f->v0->x << " " << f->v0->y << " " << f->v0->z << std::endl;
//        std::cout << f->v1->x << " " << f->v1->y << " " << f->v1->z << std::endl;
//        std::cout << f->v2->x << " " << f->v2->y << " " << f->v2->z << std::endl;

        
        Vector *v0 = f->v0;
        Vector *v1 = f->v1;
        Vector *v2 = f->v2;
        
        v0->faceList.push_back(f);
        v1->faceList.push_back(f);
        v2->faceList.push_back(f);
    }
    

    
    //Constructing the edges and the face to edge pointers (f->e = e)
    for(int i=0; i<this->faceArray.size(); i++){
        
        Face *f = &this->faceArray[i];
        
        Edge *e0 = new Edge(f->v0, f->v1);
        Edge *e1 = new Edge(f->v1, f->v2);
        Edge *e2 = new Edge(f->v2, f->v0);
        
        f->e0 = e0;
        f->e1 = e1;
        f->e2 = e2;
        
        //adding the face edges to the mesh...
        this->edgeList.push_back(*e0);
        this->edgeList.push_back(*e1);
        this->edgeList.push_back(*e2);
    }
    
    //Constructing edge to faces pointers (e->leftFace = f1 && e->rightFace = f2)
    for(int i=0; i<this->faceArray.size(); i++){
        
        Face *f = &this->faceArray[i];
        this->makeEdgeConnections(f);
    }
    
    //Finding the border edges of the mesh
    this->findOutline();

    return true;
}

bool Mesh::convertFromMshSurfToMesh(int nPts, double *pts, int nElems, int *elems){
    
    for(int i=0; i<nPts; i++){
        Vector v;
        v.x = pts[3*i+0];
        v.y = pts[3*i+1];
        v.z = pts[3*i+2];

        this->vertexArray.push_back(v);
    }
    
    for(int i=0; i<nElems; i++){
        int indV0 = elems[4*i+1];
        int indV1 = elems[4*i+2];
        int indV2 = elems[4*i+3];
        
        Vector *v0 = &this->vertexArray[indV0];
        Vector *v1 = &this->vertexArray[indV1];
        Vector *v2 = &this->vertexArray[indV2];
        
        Face f;
        f.v0 = v0;
        f.v1 = v1;
        f.v2 = v2;


        glm::vec3 c01,c02;
        Vector v01 = *f.v1 - *f.v0;
        Vector v02 = *f.v2 - *f.v0;
        Vector v01v02 = v01 % v02;
        
        c01.x = v01.x;
        c01.y = v01.y;
        c01.z = v01.z;
        
        c02.x = v02.x;
        c02.y = v02.y;
        c02.z = v02.z;
        
        
        //f.normal = glm::normalize(glm::cross(c01,c02));
        f.normal = v01v02.normalize();
        f.consistency = INT_MAX;

        this->normalList.push_back(f.center());
        //this->normalList.push_back(f.center() + glm::normalize(f.normal));
        this->normalList.push_back(f.center() + f.normal.normalize());
        this->faceArray.push_back(f);
    }

    //Constructing the backwards pointers (v->fi)
    for(int i=0; i<this->faceArray.size(); i++){
        
        Face *f = &this->faceArray[i];
        
        Vector *v0 = f->v0;
        Vector *v1 = f->v1;
        Vector *v2 = f->v2;
        
        v0->faceList.push_back(f);
        v1->faceList.push_back(f);
        v2->faceList.push_back(f);
    }

    //Constructing the edges and the face to edge pointers (f->e = e)
    for(int i=0; i<this->faceArray.size(); i++){
        
        Face *f = &this->faceArray[i];
        
        Edge *e0 = new Edge(f->v0, f->v1);
        Edge *e1 = new Edge(f->v1, f->v2);
        Edge *e2 = new Edge(f->v2, f->v0);
        
        f->e0 = e0;
        f->e1 = e1;
        f->e2 = e2;
        
        //adding the face edges to the mesh...
        this->edgeList.push_back(*e0);
        this->edgeList.push_back(*e1);
        this->edgeList.push_back(*e2);
    }
    
    //Constructing edge to faces pointers (e->leftFace = f1 && e->rightFace = f2)
    for(int i=0; i<this->faceArray.size(); i++){
        
        Face *f = &this->faceArray[i];
        this->makeEdgeConnections(f);
    }
    
    //Finding the border edges of the mesh
    this->findOutline();
    
    return true;
}

int exists(std::vector<Vector> array, Vector v){
    for(int i=0; i<array.size(); i++){
        if(v.x == array[i].x && v.y == array[i].y && v.z == array[i].z){
            return i;
        }
    }
    return -1;
}

bool Mesh::exportObj(int newVertexAmount, std::string *name){
    
    std::vector<std::string> pathVector = split(this->name, '/');
    std::string objName = pathVector[pathVector.size()-1];
    
    std::string newName;
    
    if(NULL == name) newName = "new_";
    else newName = *name;
    if(pathVector[pathVector.size()-1] != "default.obj")
        newName.append(objName);
    else
        newName = pathVector[pathVector.size()-1];
    
    std::string newPath;
    for(int i=0; i<pathVector.size()-1; i++){
        newPath.append(pathVector[i]+"/");
    }
    newPath.append(newName);
    
    std::fstream file;
    file.open(newPath, std::fstream::out);
    
    int cont = 0;
    
    for(int i=0; i<this->vertexArray.size(); i++){
        Vector *vf = &this->vertexArray[i];
        
        if(vf->toBeErased) cont++;
        else vf->id -= cont;
    }

    for(int i=0; i<this->vertexArray.size(); i++){
        if(!this->vertexArray[i].toBeErased){
            file << "v " << this->vertexArray[i].x << " " << this->vertexArray[i].y << " " << this->vertexArray[i].z << "\n";
        }
    }

    for(int i=0; i<this->faceArray.size(); i++){
        file << "f " << this->faceArray[i].v0->id + 1 << " " << this->faceArray[i].v1->id + 1<< " " << this->faceArray[i].v2->id + 1<< "\n";
    }

    file.close();
    
    if(newVertexAmount > 0){

        std::string auxPath;
        for(int i=0; i<pathVector.size()-1; i++){
            auxPath.append(pathVector[i]+"/");
        }
        auxPath.append("aux.txt");
        
        std::fstream auxFile;
        auxFile.open(auxPath, std::fstream::out);
        
        auxFile << newVertexAmount << "\n";
        
        for(int i=0; i<this->vertexArray.size(); i++){
            auxFile << this->vertexArray[i].toBeErased << "\n";
        }
        auxFile.close();
        
    }else{
        
//        std::unordered_map<int, int> indexMap;
//        int id = 1;
//        for(int i=0; i<this->vertexArray.size(); i++){
//            indexMap.insert({this->vertexArray[i].id, id++});
//        }
//
//        for(int i=0; i<this->vertexArray.size(); i++){
//            file << "v " << this->vertexArray[i].x << " " << this->vertexArray[i].y << " " << this->vertexArray[i].z << "\n";
//        }
//        
//        for(int i=0; i<this->faceArray.size(); i++){
//            file << "f "
//                 << indexMap.at(this->faceArray[i].v0->id) << " "
//                 << indexMap.at(this->faceArray[i].v1->id) << " "
//                 << indexMap.at(this->faceArray[i].v2->id) << "\n";
//        }
//        
//        file.close();
     }
    
    
    return true;
}

bool Mesh::exportHoleFacesId(std::list<int> facesArrayIndex){
    std::vector<std::string> pathVector = split(this->name, '/');
    std::string auxPath;
    for(int i=0; i<pathVector.size()-1; i++){
        auxPath.append(pathVector[i]+"/");
    }
    auxPath.append("hole.txt");
    
    std::fstream holeFile;
    holeFile.open(auxPath, std::fstream::out);
    
    for(std::list<int>::iterator i=facesArrayIndex.begin(); i!=facesArrayIndex.end(); i++){
        holeFile << *i << "\n";
    }
    
    holeFile.close();
    return true;
}

bool Mesh::getHasNormals(){
    return this->hasNormals;
}

bool Mesh::getHasTextureCoords(){
    return this->hasTextureCoords;
}

std::string Mesh::getName(){
    return this->name;
}























