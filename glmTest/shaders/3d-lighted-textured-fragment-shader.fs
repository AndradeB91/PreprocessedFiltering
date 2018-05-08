#version 150

uniform mat4 model;
uniform vec3 cameraPosition;
uniform mat3 normalMatrix;

// material settings
uniform sampler2D materialTex;
uniform float materialShininess;
uniform vec3 materialSpecularColor;

uniform sampler2D shadowMap;


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

in vec2 fragTexCoord;
in vec3 fragNormal;
in vec3 fragVert;
in vec4 ShadowCoord;

out vec4 finalColor;


vec2 poissonDisk[16] = vec2[](
                              vec2( -0.94201624, -0.39906216 ),
                              vec2( 0.94558609, -0.76890725 ),
                              vec2( -0.094184101, -0.92938870 ),
                              vec2( 0.34495938, 0.29387760 ),
                              vec2( -0.91588581, 0.45771432 ),
                              vec2( -0.81544232, -0.87912464 ),
                              vec2( -0.38277543, 0.27676845 ),
                              vec2( 0.97484398, 0.75648379 ),
                              vec2( 0.44323325, -0.97511554 ),
                              vec2( 0.53742981, -0.47373420 ),
                              vec2( -0.26496911, -0.41893023 ), 
                              vec2( 0.79197514, 0.19090188 ), 
                              vec2( -0.24188840, 0.99706507 ), 
                              vec2( -0.81409955, 0.91437590 ), 
                              vec2( 0.19984126, 0.78641367 ), 
                              vec2( 0.14383161, -0.14100790 ) 
                              );

float random(vec3 seed, int i){
    vec4 seed4 = vec4(seed,i);
    float dot_product = dot(seed4, vec4(12.9898,78.233,45.164,94.673));
    return fract(sin(dot_product) * 43758.5453);
}

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
    vec4 surfaceColor = texture(materialTex, fragTexCoord);
    vec3 surfaceToCamera = normalize(cameraPosition - surfacePos);
    
    //combine color from all the lights
    vec3 linearColor = vec3(0);
    for(int i = 0; i < numLights; ++i){
        if(i == 1)
            linearColor += ApplyLight(allLights[i], surfaceColor.rgb, normal, surfacePos, surfaceToCamera, true);
        else
            linearColor += ApplyLight(allLights[i], surfaceColor.rgb, normal, surfacePos, surfaceToCamera, false);
        
    }


//    vec3 l = normalize(allLights[0].position.xyz - surfacePos);
//    float cosTheta = clamp( dot( normal , l ), 0,1);
//    float bias = 0.00001*tan(acos(cosTheta)); // cosTheta is dot( n,l ), clamped between 0 and 1
//    bias = clamp(bias, 0,0.00003);



    
    //float bias = 0.00001; //scene 1
    float bias = 0.00003;
    
    float visibility = 1.0;
    int index;
    
//    for (int i=0;i<16;i++){
//        index = int(16.0*random(floor(fragVert.xyz*1000.0), i))%16;
//        if ( texture( shadowMap, (ShadowCoord.xy/ShadowCoord.w) + poissonDisk[index]/1000.0 ).r  <  (ShadowCoord.z-bias)/ShadowCoord.w ){
//            visibility-=0.05;
//        }
//    }
    
    for (int i=-6;i<6;i++){
        for (int j=-6;j<6;j++){
            if ( texture( shadowMap, (ShadowCoord.xy/ShadowCoord.w) + vec2(i,j)/3000.0 ).r  <  (ShadowCoord.z-bias)/ShadowCoord.w ){
                visibility-=0.005;
            }
        }
    }
    
    if(ShadowCoord.w <= 0.0f) {
        visibility = 1.0;
    } else if(ShadowCoord.x/ShadowCoord.w < 0 || ShadowCoord.y/ShadowCoord.w < 0) {
        visibility = 1.0;
    } else if(ShadowCoord.x/ShadowCoord.w >= 1 || ShadowCoord.y/ShadowCoord.w >=1) {
        visibility = 1.0;
    }
    
    //final color (after gamma correction)
    vec3 gamma = vec3(1.0/2.2);
    
    finalColor =  vec4(pow(linearColor, gamma), surfaceColor.a);
    
    finalColor.r = visibility*finalColor.r;
    finalColor.g = visibility*finalColor.g;
    finalColor.b = visibility*finalColor.b;
    
   
}


