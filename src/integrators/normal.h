/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/platform.h>
#include <core/integrator.h>

TR_NAMESPACE_BEGIN

/**
 * Surface normal integrator.
 */
struct NormalIntegrator : Integrator {
    explicit NormalIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        //v3f color(0.f, 1.f, 0.f); // Default rgb value
        SurfaceInteraction *surfaceInteraction = new SurfaceInteraction();
        if(scene.bvh->intersect(ray,*surfaceInteraction)){
            float x = surfaceInteraction->frameNs.n.x;
            float y = surfaceInteraction->frameNs.n.y;
            float z = surfaceInteraction->frameNs.n.z;
            return glm::vec3(fabs(x),fabs(y),fabs(z));
        } else{
            return glm::vec3(0.f,0.f,0.f);
        }

        // HINT: Use the scene.bvh->intersect method. It's definition is in src/accel.h
        // TODO(A1): Implement this

       //return color;
    }
};

TR_NAMESPACE_END