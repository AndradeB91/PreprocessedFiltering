//#version 150
//
//uniform mat4 camera;
//uniform mat4 model;
//
//in vec3 vert;
//in vec2 vertTexCoord;
//in vec3 vertNormal;
//
//out vec3 fragVert;
//out vec2 fragTexCoord;
//out vec3 fragNormal;
//
//void main() {
//    // Pass the tex coord straight through to the fragment shader
//    fragTexCoord = vertTexCoord;
//    fragNormal = vertNormal;
//    fragVert = vert;
//
//    // Apply all matrix transformations to the vertex
//    gl_Position = camera * model * vec4(vert, 1);
//    
//}




#version 150

uniform mat4 camera;
uniform mat4 model;
uniform mat4 DepthBiasMVP;


in vec3 vert;
in vec2 vertTexCoord;
in vec3 vertNormal;

out vec3 fragVert;
out vec2 fragTexCoord;
out vec3 fragNormal;

out vec4 ShadowCoord;

void main() {
    // Pass the tex coord straight through to the fragment shader
    fragTexCoord = vertTexCoord;
    fragNormal = vertNormal;
    fragVert = vert;
    
    // Apply all matrix transformations to the vertex
    gl_Position = camera * model * vec4(vert, 1);
    
    ShadowCoord =  DepthBiasMVP * vec4(vert, 1);
    
}