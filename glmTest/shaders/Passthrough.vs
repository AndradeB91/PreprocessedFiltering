#version 150

in vec3 vert;
in vec2 vertTexCoord;

out vec2 fragTexCoord;

void main(){

	gl_Position =  vec4(vert,1);
    fragTexCoord = (vert.xy+vec2(1,1))/2.0;
}

