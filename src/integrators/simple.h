/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Simple direct illumination integrator.
 */
struct SimpleIntegrator : Integrator {
    explicit SimpleIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
        SurfaceInteraction surfaceInteraction;
        if(scene.bvh->intersect(ray,surfaceInteraction)){
            v3f position = scene.getFirstLightPosition();
            v3f intensity = scene.getFirstLightIntensity();

            //v3f surfaceMaterial = getBSDF(surfaceInteraction)->eval(surfaceInteraction);

            surfaceInteraction.wi=position-surfaceInteraction.p;
            surfaceInteraction.wi=surfaceInteraction.frameNs.toLocal(surfaceInteraction.wi);
            surfaceInteraction.wi=glm::normalize(surfaceInteraction.wi);
            surfaceInteraction.wo=glm::normalize(surfaceInteraction.wo);

            float distance = glm::distance(surfaceInteraction.p, position);

            Ray shadowRay(surfaceInteraction.p,normalize(position-surfaceInteraction.p),Epsilon,distance);
            SurfaceInteraction shadowSurfaceInteraction;


            shadowSurfaceInteraction.wi=position-surfaceInteraction.p;
            shadowSurfaceInteraction.wi=shadowSurfaceInteraction.frameNs.toLocal(shadowSurfaceInteraction.wi);
            shadowSurfaceInteraction.wi=glm::normalize(shadowSurfaceInteraction.wi);
            shadowSurfaceInteraction.wo=glm::normalize(shadowSurfaceInteraction.wo);
            Li = (intensity/ (distance*distance)) * getBSDF(surfaceInteraction)->eval(surfaceInteraction);
            if(scene.bvh->intersect(shadowRay,shadowSurfaceInteraction)){
                Li = v3f (0.f);
            }
        }
        // TODO(A2): Implement this

        return Li;
    }
};

TR_NAMESPACE_END