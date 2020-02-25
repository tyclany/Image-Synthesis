/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator with MIS
 */
struct DirectIntegrator : Integrator {
    explicit DirectIntegrator(const Scene& scene) : Integrator(scene) {
        m_emitterSamples = scene.config.integratorSettings.di.emitterSamples;
        m_bsdfSamples = scene.config.integratorSettings.di.bsdfSamples;
        m_samplingStrategy = scene.config.integratorSettings.di.samplingStrategy;
    }

    static inline float balanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
        float f = nf * fPdf, g = ng * gPdf;
        return f / (f + g);
    }

    void sampleSphereByCosineHemisphere(const p2f& sample,
                                        const v3f& n,
                                        const p3f& pShading,
                                        const v3f& emitterCenter,
                                        float emitterRadius,
                                        v3f& wiW,
                                        float& pdf) const {
        wiW = Warp::squareToCosineHemisphere(sample);
        pdf = Warp::squareToCosineHemispherePdf(wiW);
        // TODO(A3): Implement this
    }

    void sampleSphereByArea(const p2f& sample,
                            const p3f& pShading,
                            const v3f& emitterCenter,
                            float emitterRadius,
                            v3f& pos,
                            v3f& ne,
                            v3f& wiW,
                            float& pdf) const {
        ne = Warp::squareToUniformSphere(sample);
        v3f samplePoint = ne * emitterRadius;
        pos = samplePoint + emitterCenter;
        wiW = glm::normalize(pos - pShading);
        pdf = INV_FOURPI / pow(emitterRadius, 2);
        // TODO(A3): Implement this
    }

    void sampleSphereBySolidAngle(const p2f& sample,
                                  const p3f& pShading,
                                  const v3f& emitterCenter,
                                  float emitterRadius,
                                  v3f& wiW,
                                  float& pdf) const {
        v3f wc = normalize(emitterCenter - pShading);
        Frame sphereFrame = Frame(wc);
        float sinThetaMax2 = emitterRadius * emitterRadius / glm::pow(glm::distance(pShading, emitterCenter), 2);
        float cosThetaMax = sqrt(glm::max(0.f, 1 - sinThetaMax2));
        wiW = Warp::squareToUniformCone(sample, cosThetaMax);
        wiW = glm::normalize(sphereFrame.toWorld(wiW));
        pdf = Warp::squareToUniformConePdf(cosThetaMax);

        // TODO(A3): Implement this
    }

    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        SurfaceInteraction surfaceInteraction = SurfaceInteraction();
        SurfaceInteraction shadowSurfaceInteraction = SurfaceInteraction();
        if (scene.bvh->intersect(ray, surfaceInteraction)){
            v3f emission = getEmission(surfaceInteraction);
            if (emission != v3f(0.f)) {
                return emission;
            }
            for(int i = 0; i<m_emitterSamples;i++){
                float emPdf;
                size_t id = selectEmitter(sampler.next(), emPdf);
                const Emitter& em = getEmitterByID(id);
                v3f emCenter = scene.getShapeCenter(em.shapeID);
                float emRadius = scene.getShapeRadius(em.shapeID);
                float pdf;
                p3f& pShading = surfaceInteraction.p;
                v3f pos;
                v3f ne;
                v3f wiW;
                sampleSphereByArea(sampler.next2D(),pShading,emCenter,emRadius,pos,ne,wiW,pdf);
                float cosThetaO = glm::dot(ne, wiW);
                float distance = glm::distance(pShading, pos);
                Ray shadowRay = Ray(pShading, wiW);
                if (scene.bvh->intersect(shadowRay, shadowSurfaceInteraction)) {
                    v3f emission = getEmission(shadowSurfaceInteraction);
                    if (pdf != 0) {
                        v3f G(0.f);
                        if(glm::dot(pos-pShading,ne)>=0){
                            G = emission * cosThetaO / glm::pow(distance, 2);
                        }
                        surfaceInteraction.wi = surfaceInteraction.frameNs.toLocal(wiW);
                        v3f brdf = getBSDF(surfaceInteraction)->eval(surfaceInteraction);
                        Lr += brdf * G / pdf / emPdf;
                    }
                }
            }
        }

        // TODO(A3): Implement this

        return Lr/m_emitterSamples;
    }

    v3f renderCosineHemisphere(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        SurfaceInteraction surfaceInteraction = SurfaceInteraction();
        SurfaceInteraction shadowSurfaceInteraction = SurfaceInteraction();
        float emPdf;
        size_t id = selectEmitter(sampler.next(), emPdf);
        const Emitter& em = getEmitterByID(id);
        v3f emCenter = scene.getShapeCenter(em.shapeID);
        float emRadius = scene.getShapeRadius(em.shapeID);
        float pdf;
        p3f pShading;
        v3f wiW;
        if (scene.bvh->intersect(ray, surfaceInteraction)){
            v3f emission = getEmission(surfaceInteraction);
            if (emission != v3f(0.f)) {
                return emission;
            }
            for(int i = 0; i<m_bsdfSamples;i++){
                sampleSphereByCosineHemisphere(sampler.next2D(), surfaceInteraction.p, pShading, emCenter, emRadius, wiW, pdf);
                surfaceInteraction.wi=wiW;
                wiW = glm::normalize(surfaceInteraction.frameNs.toWorld(wiW));
                Ray shadowRay = Ray(surfaceInteraction.p, wiW);
                if (scene.bvh->intersect(shadowRay, shadowSurfaceInteraction)) {
                    v3f emission = getEmission(shadowSurfaceInteraction);
                    if (pdf != 0.f) {
                        v3f brdf = getBSDF(surfaceInteraction)->eval(surfaceInteraction);
                        Lr += brdf * emission / pdf;
                    }
                }
            }
        }
        // TODO(A3): Implement this
        return Lr/m_bsdfSamples;
    }

    v3f renderBSDF(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        SurfaceInteraction surfaceInteraction;
        SurfaceInteraction shadowSurfaceInteraction;
        if (scene.bvh->intersect(ray, surfaceInteraction)) {
            v3f emission = getEmission(surfaceInteraction);
            if (emission != v3f(0.f)) {
                return emission;
            }
            for (int j = 0; j <m_bsdfSamples; j++) {
                float pdf;
                v3f brdf = getBSDF(surfaceInteraction)->sample(surfaceInteraction, sampler, &pdf);
                v3f wi = surfaceInteraction.wi;
                wi = surfaceInteraction.frameNs.toWorld(wi);
                Ray shadowRay = Ray(surfaceInteraction.p, wi);
                if (scene.bvh->intersect(shadowRay, shadowSurfaceInteraction)) {
                    v3f emission = getEmission(shadowSurfaceInteraction);
                    if(pdf!=0){
                        Lr += brdf * emission/pdf;
                    }
                }
            }
        }
        // TODO(A3): Implement this

        return Lr/m_bsdfSamples;
    }

    v3f renderSolidAngle(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        SurfaceInteraction surfaceInteraction;
        SurfaceInteraction shadowSurfaceInteraction;
        if(scene.bvh->intersect(ray, surfaceInteraction)){
            v3f emission = getEmission(surfaceInteraction);
            if (emission != v3f(0.f)) {
                return emission;
            }

            for(int i = 0; i < m_emitterSamples; i++){
                float emPdf;
                size_t id = selectEmitter(sampler.next(), emPdf);
                const Emitter& em = getEmitterByID(id);
                v3f emCenter = scene.getShapeCenter(em.shapeID);
                float emRadius = scene.getShapeRadius(em.shapeID);
                float pdf;
                v3f wiW;
                sampleSphereBySolidAngle(sampler.next2D(), surfaceInteraction.p, emCenter, emRadius, wiW, pdf);
                Ray shadow_ray = Ray(surfaceInteraction.p, normalize(wiW));
                if(scene.bvh->intersect(shadow_ray, shadowSurfaceInteraction)){
                    v3f emission = getEmission(shadowSurfaceInteraction);
                    if(pdf!=0){
                        surfaceInteraction.wi = normalize(surfaceInteraction.frameNs.toLocal(wiW));
                        v3f brdf = getBSDF(surfaceInteraction)->eval(surfaceInteraction);
                        Lr += brdf * emission / pdf / emPdf;
                    }
                }
            }
        }

        return Lr/m_emitterSamples;
        // TODO(A3): Implement this
    }

    v3f renderMIS(const Ray& ray, Sampler& sampler) const {

        v3f Lr(0.f);
        v3f LrEmitter(0.f);
        v3f LrBSDF(0.f);

        SurfaceInteraction surfaceInteraction;
        SurfaceInteraction shadowSurfaceInteraction;
        if(scene.bvh->intersect(ray, surfaceInteraction)){
            v3f emission = getEmission(surfaceInteraction);
            if (emission != v3f(0.f)) {
                return emission;
            }

            for(int i = 0; i < m_emitterSamples; i++){
                float emPdf;
                size_t id = selectEmitter(sampler.next(), emPdf);
                const Emitter& em = getEmitterByID(id);
                v3f emCenter = scene.getShapeCenter(em.shapeID);
                float emRadius = scene.getShapeRadius(em.shapeID);
                float pdf;
                float BSDFpdf;
                v3f wiW;
                sampleSphereBySolidAngle(sampler.next2D(), surfaceInteraction.p, emCenter, emRadius, wiW, pdf);
                v3f BSDFbrdf = getBSDF(surfaceInteraction)->sample(surfaceInteraction, sampler, &BSDFpdf);
                Ray shadow_ray = Ray(surfaceInteraction.p, normalize(wiW));
                if(scene.bvh->intersect(shadow_ray, shadowSurfaceInteraction)){
                    v3f emission = getEmission(shadowSurfaceInteraction);
                    if(pdf!=0){
                        surfaceInteraction.wi = normalize(surfaceInteraction.frameNs.toLocal(wiW));
                        float weightForDiffuse = balanceHeuristic(m_emitterSamples,pdf*emPdf,m_bsdfSamples,BSDFpdf);
                        v3f brdf = getBSDF(surfaceInteraction)->eval(surfaceInteraction);
                        LrEmitter += brdf * emission *weightForDiffuse / pdf / emPdf;
                    }
                }
            }

            for (int j = 0; j <m_bsdfSamples; j++) {
                float pdf;
                v3f brdf = getBSDF(surfaceInteraction)->sample(surfaceInteraction, sampler, &pdf);
                v3f wi = surfaceInteraction.wi;
                wi = surfaceInteraction.frameNs.toWorld(wi);
                Ray shadowRay = Ray(surfaceInteraction.p, wi);
                float emPdf;
                size_t id = selectEmitter(sampler.next(), emPdf);
                const Emitter& em = getEmitterByID(id);
                v3f emCenter = scene.getShapeCenter(em.shapeID);
                float emRadius = scene.getShapeRadius(em.shapeID);
                float pfpdf;
                v3f wiW;
                sampleSphereBySolidAngle(sampler.next2D(), surfaceInteraction.p, emCenter, emRadius, wiW, pfpdf);
                if (scene.bvh->intersect(shadowRay, shadowSurfaceInteraction)) {
                    v3f emission = getEmission(shadowSurfaceInteraction);
                    if(pdf!=0&&pfpdf!=0){
                        float weightForpf=balanceHeuristic(m_bsdfSamples, pdf, m_emitterSamples, pfpdf);
                        LrBSDF += brdf * emission*weightForpf/pdf;
                    }
                }
            }
        }
        if (m_emitterSamples == 0) {
            Lr = LrBSDF / m_bsdfSamples;
        }
        else if (m_bsdfSamples == 0) {
            Lr = LrEmitter / m_emitterSamples;
        }
        else {
            Lr = LrEmitter / m_emitterSamples + LrBSDF / m_bsdfSamples;
        }

        return Lr;
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        if (m_samplingStrategy == ESamplingStrategy::EMIS)
            return this->renderMIS(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::EArea)
            return this->renderArea(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::ESolidAngle)
            return this->renderSolidAngle(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::ECosineHemisphere)
            return this->renderCosineHemisphere(ray, sampler);
        else
            return this->renderBSDF(ray, sampler);
    }

    size_t m_emitterSamples;     // Number of emitter samples
    size_t m_bsdfSamples;        // Number of BSDF samples
    ESamplingStrategy m_samplingStrategy;   // Sampling strategy to use
};

TR_NAMESPACE_END