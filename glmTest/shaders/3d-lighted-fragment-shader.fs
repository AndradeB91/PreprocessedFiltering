#version 150

uniform mat4 model;
uniform vec3 cameraPosition;
uniform mat3 normalMatrix;

// material settings
uniform float materialShininess;
uniform vec3 materialSpecularColor;

#define MAX_LIGHTS 10
uniform int numLights;

uniform struct Light {
    vec4 position;
    vec3 intensities; //a.k.a the color of the light
    float attenuation;
    float ambientCoefficient;
    float coneAngle;
    vec3 coneDirection;
} allLights[MAX_LIGHTS];

in vec3 fragNormal;
in vec3 fragVert;

out vec4 finalColor;


vec3 ApplyLight(Light light, vec3 surfaceColor, vec3 normal, vec3 surfacePos, vec3 surfaceToCamera, bool flashLight) {
    vec3 surfaceToLight;
    float attenuation = 1.0;
    
    float attenuationAdd = 0.0;
    if(light.position.w == 0.0) {
        //directional light
        surfaceToLight = normalize(light.position.xyz - surfacePos);
        attenuation = 1.0; //no attenuation for directional lights
    } else {
        //point light
        surfaceToLight = normalize(light.position.xyz - surfacePos);
        float distanceToLight = length(light.position.xyz - surfacePos);
        attenuation = 1.0 / (1.0 + light.attenuation * pow(distanceToLight, 2));

        //cone restrictions (affects attenuation)
        float lightToSurfaceAngle = degrees(acos(dot(-surfaceToLight, normalize(light.coneDirection))));
        if(lightToSurfaceAngle < 10){
            attenuationAdd = (11-lightToSurfaceAngle);
        }
        if(lightToSurfaceAngle > light.coneAngle){
            float difAngles = lightToSurfaceAngle - light.coneAngle;
            
            if(flashLight){
                attenuation /= (1+difAngles*3);

            }
            else{
                attenuation /= (1+difAngles*5);
            }

        }
    }

    //ambient
    vec3 ambient = light.ambientCoefficient * surfaceColor.rgb * light.intensities;

    //diffuse
    float diffuseCoefficient = max(0.0, dot(normal, surfaceToLight));
    vec3 diffuse = diffuseCoefficient * surfaceColor.rgb * light.intensities;

    //specular
    float specularCoefficient = 0.0;
    if(diffuseCoefficient > 0.0)
        specularCoefficient = pow(max(0.0, dot(surfaceToCamera, reflect(-surfaceToLight, normal))), materialShininess);
    vec3 specular = specularCoefficient * materialSpecularColor * light.intensities;

    //linear color (color before gamma correction)
    return ambient + attenuation*(diffuse + specular);
}

void main(){
    vec3 normal = normalize(normalMatrix * fragNormal);
    vec3 surfacePos = vec3(model * vec4(fragVert, 1));
    vec4 surfaceColor = vec4(1,1,1,1);
    vec3 surfaceToCamera = normalize(cameraPosition - surfacePos);

    //combine color from all the lights
    vec3 linearColor = vec3(0);
    for(int i = 0; i < numLights; ++i){
        if(i == 1)
            linearColor += ApplyLight(allLights[i], surfaceColor.rgb, normal, surfacePos, surfaceToCamera, true);
        else
            linearColor += ApplyLight(allLights[i], surfaceColor.rgb, normal, surfacePos, surfaceToCamera, false);

    }

    //final color (after gamma correction)
    vec3 gamma = vec3(1.0/2.2);
    finalColor = vec4(pow(linearColor, gamma), surfaceColor.a);
}













