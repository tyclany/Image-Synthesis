
#version 330 core

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 normalMat;


layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;

out vec3 vNormal;
out vec3 vPos;
void main()
{
    gl_Position = projection*view*model*vec4(position,1.0);
    vPos = vec3(model*vec4(position,1.0));
    vNormal = vec3(normalMat*vec4(normal,0)); // cast the result to vec3

}