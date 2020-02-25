/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once
#include <random>
#include "bsdfs/phong.h"

TR_NAMESPACE_BEGIN

/**
 * Reflection occlusion integrator
 */
struct ROIntegrator : Integrator {

    float m_exponent;

    explicit ROIntegrator(const Scene& scene) : Integrator(scene) {
        m_exponent = scene.config.integratorSettings.ro.exponent;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }


    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
        SurfaceInteraction surfaceInteraction;
        SurfaceInteraction shadowSurfaceInteraction;
        if (scene.bvh->intersect(ray, surfaceInteraction)){
            v3f x = surfaceInteraction.p;
            v3f	wi = Warp::squareToPhongLobe(sampler.next2D(), m_exponent);
            float pdf = Warp::squareToPhongLobePdf(wi, m_exponent);
            float cosAlpha = wi.z;
            v3f wo = surfaceInteraction.wo;
            v3f wr = reflect(wo);

            Frame reflectFrame = Frame(surfaceInteraction.frameNs.toWorld(wr));
            wi = reflectFrame.toWorld(wi);
            wi = surfaceInteraction.frameNs.toLocal(wi);
            float cosTheta = fmax(wi.z, 0);
            wi = surfaceInteraction.frameNs.toWorld(wi);
            Ray shadowRay = Ray(x, wi);
            if (scene.bvh->intersect(shadowRay, shadowSurfaceInteraction)) {
                return Li;
            }
            else {
                Li = v3f((m_exponent + 2) *INV_TWOPI*fmax(pow(cosAlpha, m_exponent), 0)*cosTheta / pdf);
                return Li;
            }

        }
    }
        // TODO(A3): Implement this
};

TR_NAMESPACE_END