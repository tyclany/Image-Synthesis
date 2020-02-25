


#version 330 core

#define M_PI 3.14159265358979323846f
#define INV_TWOPI  0.15915494309189533577f


uniform vec3 camPos;
uniform vec3 lightPos;
uniform vec3 lightIntensity;
uniform vec3 rho_d;
uniform vec3 rho_s;
uniform float exponent;

in vec3 vNormal;
in vec3 vPos;

out vec3 phongFragmentColor;


void main()
{
      vec3 diffuseColor = rho_d/M_PI * lightIntensity;
      float distance = distance(lightPos,vPos);
      float rvProduct= dot(normalize(reflect((lightPos-vPos),vNormal)),normalize(camPos-vPos));
      vec3 specularColor = rho_s*lightIntensity*(exponent+2)*INV_TWOPI * max(0.f,pow(rvProduct,exponent));
      phongFragmentColor = (diffuseColor+specularColor)/(distance*distance) * max(0.f, dot(vNormal,normalize(lightPos-vPos)));
}