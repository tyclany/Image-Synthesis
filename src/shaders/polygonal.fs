/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core


#define PI 3.14159265359
#define MAX_NUM_EMITTER_TRIANGLES 40 // Max number of emitter triangles allowed (tuned for A4)
uniform float emitterVertices[MAX_NUM_EMITTER_TRIANGLES*3*3]; // Need to specify a size at compile time (max size is 512 floats)

uniform int nbTriangles; // Use nbTriangles to index the buffer correctly
uniform vec3 lightIrradiance;
uniform vec3 albedo;
uniform vec2 windowSize; // [width, height]

uniform sampler2D cvTerm; // Creates a 2D texture we can sample from (range x,y = [0,1])

in vec3 vNormal;
in vec3 vPos;

out vec3 color;

// Compute edge (v1--v2) contribution
float getEdgeContrib(vec3 v1, vec3 v2, vec3 pos) {
	// Adapt your getEdgeContrib code from the offline part
	float value = 0.f;
	float theta = acos(dot(normalize(v1-pos),normalize(v2-pos)));
	vec3 yi = normalize(cross(normalize(v2-pos),normalize(v1-pos)));
	value = theta*(dot(yi,vNormal));
// TODO(A4): Implement this
	return value;
}


void main()
{	
	// 1) Extract vertices of triangles from `emitterVertices` buffer using `nbTriangles`
	// 2) Calculate G term
	// 3) Subtract modification term for G after extracting it from texture (use built-in `texture()` function)
	//	    e.g. `vec3 delta = texture(cvTerm, coords).xyz;`
    float sum = 0.f;
    for(int i = 0 ; i < nbTriangles; i++){
        vec3 v1 = vec3 (emitterVertices[9*i],emitterVertices[9*i+1],emitterVertices[9*i+2]);
        vec3 v2 = vec3 (emitterVertices[9*i+3],emitterVertices[9*i+4],emitterVertices[9*i+5]);
        vec3 v3 = vec3 (emitterVertices[9*i+6],emitterVertices[9*i+7],emitterVertices[9*i+8]);
        float contrib0 = getEdgeContrib(v1,v2,vPos);
        float contrib1 = getEdgeContrib(v2,v3,vPos);
        float contrib2 = getEdgeContrib(v3,v1,vPos);
        sum = sum+contrib0+contrib1+contrib2;
    }
    vec3 G = albedo * lightIrradiance * sum / (2*PI * PI);
    float scaledX = gl_FragCoord.x/(windowSize[0]);
    float scaledY = gl_FragCoord.y/(windowSize[1]);
    vec2 scaledPosition = vec2(scaledX,scaledY);
    vec3 delta = texture(cvTerm, scaledPosition).xyz;
	color = vec3(1.0, 0 , 0);
	color = G-delta;

    // TODO(A4): Implement this
}

