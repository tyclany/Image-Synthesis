/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Path tracer integrator
 */
struct PathTracerIntegrator : Integrator {
    explicit PathTracerIntegrator(const Scene& scene) : Integrator(scene) {
        m_isExplicit = scene.config.integratorSettings.pt.isExplicit;
        m_maxDepth = scene.config.integratorSettings.pt.maxDepth;
        m_rrDepth = scene.config.integratorSettings.pt.rrDepth;
        m_rrProb = scene.config.integratorSettings.pt.rrProb;
    }


    v3f renderImplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
        v3f Li(0.f);
        v3f acc = v3f(1.f);
        int depth = m_maxDepth;
        v3f emission = getEmission(hit);

        if(emission != v3f(0.f)&& glm::dot(v3f(0, 0, 1), hit.wo)> 0){
            return acc*emission;
        }
        while (depth != -1) {
            emission = getEmission(hit);
            if (emission != v3f(0.f) && glm::dot(v3f(0, 0, 1), hit.wo) > 0) {
                acc *= emission;
                return acc;
            } else {
                float pdf;
                v3f bsdf = getBSDF(hit)->sample(hit, sampler, &pdf);
                acc *= bsdf/pdf;
                v3f wi = glm::normalize(hit.frameNs.toWorld(hit.wi));
                Ray next = Ray(hit.p, wi);
                if (scene.bvh->intersect(next, hit)) {
                    depth--;
                } else {
                    return v3f(0.f);
                }
            }
        }
        // TODO(A5): Implement this
        return Li;
    }
    v3f renderExplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& duplicate) const {
        v3f Li(0.f);
        v3f brdfAcc = v3f(1.f);
        v3f accDir = v3f(0.f);
        int depth = m_maxDepth;
        SurfaceInteraction hit = duplicate;
        v3f emission = getEmission(hit);
        if(emission != v3f(0.f)){
            return emission;
        }

        while(depth!=0){
            float emPdf;
            size_t id = selectEmitter(sampler.next(), emPdf);
            const Emitter& em = getEmitterByID(id);
            v3f n, pos;
            float pdf;
            sampleEmitterPosition(sampler, em, n, pos, pdf);
            v3f wi = pos - hit.p;
            wi = glm::normalize(wi);

            hit.wi = hit.frameNs.toLocal(wi);
            float cosTheta = glm::dot(-wi, n);

            if (cosTheta < 0) {
                cosTheta = 0;
            }

            Ray shadowRay = Ray(hit.p, wi);
            SurfaceInteraction shadowSurfaceInteraction;
            if (scene.bvh->intersect(shadowRay, shadowSurfaceInteraction)) {
                v3f bsdf = getBSDF(hit)->eval(hit);
                v3f shadowEmission = getEmission(shadowSurfaceInteraction);
                if (shadowEmission != v3f(0.f)) {
                    float G = cosTheta / (glm::distance2(pos, hit.p));
                    accDir += brdfAcc * bsdf * shadowEmission * G / pdf / emPdf;
                }
            }
            float brdfPdf = 0.f;
            v3f brdf = v3f(0.f);
            SurfaceInteraction buffer = hit;
            do{
                brdf = getBSDF(hit)->sample(hit, sampler, &brdfPdf);
                v3f direction = hit.wi;
                direction = hit.frameNs.toWorld(direction);
                direction = glm::normalize(direction);
                Ray nextRay(hit.p, direction);
                bool hitSurface = scene.bvh->intersect(nextRay, buffer);
                if (!hitSurface) {
                    return accDir;
                }
            }while(getEmission(buffer)!=v3f(0.f));
            if(brdfPdf!=0){
                brdfAcc *= brdf/brdfPdf;
            }
            hit = buffer;
            if(depth <= m_maxDepth-m_rrDepth){
                if(sampler.next()>m_rrProb){
                    return accDir;
                }else{
                    brdfAcc /= m_rrProb;
                }
            }
            depth --;
        }
        return accDir;
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        Ray r = ray;
        SurfaceInteraction hit;

        if (scene.bvh->intersect(r, hit)) {
            if (m_isExplicit)
                return this->renderExplicit(ray, sampler, hit);
            else
                return this->renderImplicit(ray, sampler, hit);
        }
        return v3f(0.0);
    }

    int m_maxDepth;     // Maximum number of bounces
    int m_rrDepth;      // When to start Russian roulette
    float m_rrProb;     // Russian roulette probability
    bool m_isExplicit;  // Implicit or explicit
};

TR_NAMESPACE_END
