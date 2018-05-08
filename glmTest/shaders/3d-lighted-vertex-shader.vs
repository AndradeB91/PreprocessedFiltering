#version 150

uniform mat4 camera;
uniform mat4 model;

in vec3 vert;

out vec3 fragVert;

void main() {
    // Pass the tex coord straight through to the fragment shader
    fragNormal = vertNormal;
    fragVert = vert;

    // Apply all matrix transformations to the vertex
    gl_Position = camera * model * vec4(vert, 1);
}