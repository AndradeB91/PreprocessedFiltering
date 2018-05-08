#version 150

uniform mat4 depthMVP;

in vec3 vert;
in vec3 vertNormal;

out vec3 fragNormal;

void main() {
    
    fragNormal = vertNormal; //ESFERAS APARECIAM PQ N TINHA ESSA LINHA AQUI
    gl_Position = depthMVP * vec4(vert, 1);

}