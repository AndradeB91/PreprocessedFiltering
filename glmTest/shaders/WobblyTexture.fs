#version 150

in vec2 fragTexCoord;

uniform sampler2D texture2D;

out vec4 finalColor;

void main(){

    vec4 color = texture ( texture2D, fragTexCoord);
    float z = texture( texture2D, fragTexCoord ).r;

    float n = 0.2;
    float f = 10000.0;
    float c = (2.0 * n) / (f + n - z * (f - n));
    finalColor.rgb = vec3(c);


    //finalColor.rgb = color.rgb;
}


