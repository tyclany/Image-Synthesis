


#version 330 core
#define M_PI 3.14159265358979323846f


uniform vec3 camPos;
uniform vec3 lightPos;
uniform vec3 lightIntensity;
uniform vec3 albedo;

in vec3 vNormal;
in vec3 vPos;
out vec3 diffuseColor;

void main()
{
    float productNormalDistance = dot(vNormal,normalize(lightPos-vPos));
    float distance = distance(lightPos,vPos);
    diffuseColor = albedo/M_PI * lightIntensity/(distance*distance) * max(0.f, productNormalDistance);

}