/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Ambient occlusion integrator
 */
struct AOIntegrator : Integrator {

	// Use this in your switch statement to select the sampling type 
	ESamplingType m_samplingStrategy;

    explicit AOIntegrator(const Scene& scene) : Integrator(scene) { 
		m_samplingStrategy = scene.config.integratorSettings.ao.sampling_type;
	}

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
        SurfaceInteraction surfaceInteraction;
        SurfaceInteraction shadowSurfaceInteraction;

        if(scene.bvh->intersect(ray,surfaceInteraction)){
            v3f x = surfaceInteraction.p;
            v3f wi;
            float pdf;

            if (m_samplingStrategy == ESpherical) {
                wi = Warp::squareToUniformSphere(sampler.next2D());
                pdf = Warp::squareToUniformSpherePdf();
            }
            else if (m_samplingStrategy == EHemispherical ) {
                wi = Warp::squareToUniformHemisphere(sampler.next2D());
                pdf = Warp::squareToUniformHemispherePdf(wi);
            }
            else if (m_samplingStrategy == ECosineHemispherical){
                wi = Warp::squareToCosineHemisphere(sampler.next2D());
                pdf = Warp::squareToCosineHemispherePdf(wi);
            }
            float cosTheta = glm::max(wi.z, 0.f);
            wi = glm::normalize(surfaceInteraction.frameNs.toWorld(wi));
            Ray shadowRay = Ray(x, wi);
            shadowRay.max_t = scene.aabb.getBSphere().radius/2;
            if (scene.bvh->intersect(shadowRay, shadowSurfaceInteraction)) {
                return Li;
            }
            else {
                Li = v3f(INV_PI*cosTheta / pdf);
                return Li;
            }
        }

		/*
		Use the m_sampling_type variable to set wi and the corresponding pdf 
		appropriately for sphere, hemisphere, or cosine sampling.

		You can use a switch statement or an if/else block.

		The m_sampling_type variable is an enum. The different values of the enum 
		can be accessed through:
		ESamplingType::ESpherical
		ESamplingType::EHemispherical
		ESamplingType::ECosineHemispherical
		*/
		
        // TODO(A3): Implement this
    }
};

TR_NAMESPACE_END