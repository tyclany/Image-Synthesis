/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <tiny_obj_loader.h>
#define RAY_EPS_CV 1e-5 // Use when setting min and max dist for ray in control variates code
TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator for polygonal light sources
 * Follows Arvo '94.
 */
struct PolygonalIntegrator : Integrator {

	float m_alpha;             // Control variates "strength"
	size_t m_visSamples;       // # of samples to estimate h - alpha*g
	bool m_traceShadows;       // Trace shadows or not
	EPolygonalMethod m_method; // Method to use (Arvo, or control variates)

	std::vector<std::vector<v3f>> m_triangles; // Data structure to store triangles

    explicit PolygonalIntegrator(const Scene& scene) : Integrator(scene) {
        m_alpha = scene.config.integratorSettings.poly.alpha;
        m_visSamples = scene.config.integratorSettings.poly.visSamples;
        m_traceShadows = scene.config.integratorSettings.poly.traceShadows;
        m_method = scene.config.integratorSettings.poly.method;

		/**
		 * 1) Get # of triangles on emitter
		 * 2) Store vertices in m_triangles
		 */
// TODO(A4): Implement this
        size_t shapeID = scene.getFirstLight();
        auto shape = scene.worldData.shapes[shapeID];
        int size = shape.mesh.indices.size();
        for(int i =0;i<size;i=i+3){
            //v3f position = scene.getObjectVertexPosition(shapeID,0);
            std::vector<v3f> tri(3);    // Create empty triangle (3 vertices / triangle)
            tri[0] = scene.getObjectVertexPosition(shapeID,i);      // Add vertices to triangle
            tri[1] = scene.getObjectVertexPosition(shapeID,i+1);
            tri[2] = scene.getObjectVertexPosition(shapeID,i+2);
            m_triangles.push_back(tri); // Append triangle vector to vector
        }
    }

    /// Reflect
    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    /**
     * === PHONG BONUS ONLY ===
     * Compute the following integral:
     *    T(a, b, n, x) = \int_0^x [a \cos(\theta) + b \sin(\theta)]ˆn d\theta
     * Uses a recurrent relation (see Snyder's note, 1996)
     *
     * Series function:
     *    T_sum(a, b, n, x) = \sum_{i=0}ˆ{(n-1)/2} T(a, b, 2i+1, x)
     * assuming n is _odd_
     */
    float cosineSinePowerIntegralSum(float a, float b, int exp, float theta) const {
        if (exp % 2 == 0) exp += 1; // Make exponent odd
        float Tsum = 0.f;

		// Implementing this function may be useful if you attempt the bonus

        // TODO(A4): Implement this

        return Tsum;
    }

    /**
     * Compute edge (v1--v2) contribution
	 * The exp term is only needed if you attempt the bonus, otherwise, you can ignore it
     */
    float getEdgeContrib(const v3f& v1, const v3f& v2, const SurfaceInteraction& i, int exp = 0) const {
        float contrib = 0.f;
        float theta = acos(clamp(glm::dot( glm::normalize(v1-i.p),glm::normalize(v2-i.p)),-1.f,1.f));
        v3f yi = glm::normalize(cross(glm::normalize(v2-i.p),glm::normalize(v1-i.p)));
        contrib = theta*(glm::dot(yi,i.frameNs.n));
        // TODO(A4): Implement this

        return contrib;
    }
	   

    /// Direct illumination using Arvo '94 analytic solution for polygonal lights
    v3f renderAnalytic(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        v3f rho(0.f);
        float sum =0.f ;
        SurfaceInteraction surfaceInteraction;
        if(scene.bvh->intersect(ray, surfaceInteraction)) {
            v3f emission = getEmission(surfaceInteraction);
            if (emission != v3f(0.f)) {
                return emission;
            }
            for(int i = 0; i<m_triangles.size() ;i++){
                float contrib0 = getEdgeContrib(m_triangles[i][0],m_triangles[i][1],surfaceInteraction);
                float contrib1 = getEdgeContrib(m_triangles[i][1],m_triangles[i][2],surfaceInteraction);
                float contrib2 = getEdgeContrib(m_triangles[i][2],m_triangles[i][0],surfaceInteraction);
                sum += contrib0+contrib1+contrib2;
            }

            size_t shapeID = scene.getFirstLight();
            const Emitter& em = getEmitterByID(getEmitterIDByShapeID(shapeID));
            surfaceInteraction.wi = v3f(0,0,1);
            rho = getBSDF(surfaceInteraction)->eval(surfaceInteraction);
            Lr += rho * em.getPower() * sum * INV_TWOPI /em.area;
        }
            // TODO(A4): Implement this

        return Lr;
    }

    /**
     * Stand-alone estimator for h - alpha*g (with primary ray)
     * Trace a primary ray, check for emitter hit, and then call `estimateVisDiff()`
     * Used by polygonal render pass
     */
    v3f estimateVisDiffRealTime(const Ray& ray, Sampler& sampler, const Emitter& em) {
        v3f D(0.f);

        SurfaceInteraction hit;
        if (!scene.bvh->intersect(ray, hit)) return D;

        const BSDF* bsdf = getBSDF(hit);
        if (bsdf->isEmissive()) return D;

        hit.wi = v3f(0, 0, 1); // Trick to get 1/pi * albedo without cosine term
        D = estimateVisDiff(sampler, hit, em);

        return D;
    }

    /// Stand-alone estimator for h - alpha*g (without primary ray)
	/// Use RAY_EPS_CV when setting min and max dist for shadow ray
    v3f estimateVisDiff(Sampler& sampler, SurfaceInteraction& i, const Emitter& em) const {
        v3f secondI=v3f(0.f);
        for(int j=0; j<m_visSamples;j++){
            v3f h_term = v3f(0.f);
            v3f g_term = v3f(0.f);
            v3f pe;    // Point on emitter
            v3f ne;    // Surface normal at point
            float pdf; // PDF of choosing point
            sampleEmitterPosition(sampler, em, ne, pe, pdf);
            SurfaceInteraction shadowSurfaceInteraction;
            v3f wiW = glm::normalize(pe-i.p);
            i.wi = i.frameNs.toLocal(wiW);
            float cosThetaO = glm::dot(ne, -wiW);
            float distance = glm::distance(i.p, pe);
            Ray shadowRay = Ray(i.p, wiW,RAY_EPS_CV);
            double G(0.f);
            v3f brdf = getBSDF(i)->eval(i);
            G = glm::max(cosThetaO,0.f) / glm::pow(distance, 2);
            if(scene.bvh->intersect(shadowRay, shadowSurfaceInteraction)){
                v3f emission = getEmission(shadowSurfaceInteraction);
                h_term = brdf* emission*G/pdf;
            }
            g_term = brdf* G*em.getRadiance()/pdf;
            secondI += h_term - m_alpha * g_term;
        }
        secondI = secondI/m_visSamples;
        // TODO(A4): Implement this
        return secondI;
    }

    /// Control variates using Arvo '94 for direct illumination; ray trace shadows

    v3f renderControlVariates(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        v3f secondI(0.f);
        v3f G(0.f);
        SurfaceInteraction surfaceInteraction;
        //const Emitter& em = getEmitterByID(getEmitterIDByShapeID(scene.getFirstLight()));
        // TODO(A4): Implement this
        if(scene.bvh->intersect(ray,surfaceInteraction)){
            v3f emission = getEmission(surfaceInteraction);
            if (emission != v3f(0.f)) {
                return emission;
            }
            G=renderAnalytic(ray,sampler);
            secondI = estimateVisDiff(sampler,surfaceInteraction,getEmitterByID(getEmitterIDByShapeID(scene.getFirstLight())));
        }
        Lr = m_alpha*G+secondI;
        Lr = clampBelow(Lr,0.f);
        return Lr;
    }

    /// Direct illumination using surface area sampling
    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        SurfaceInteraction surfaceInteraction;
        SurfaceInteraction shadowSurfaceInteracion;
        size_t shapeID = scene.getFirstLight();
        const Emitter& em = getEmitterByID(getEmitterIDByShapeID(shapeID));
        v3f Lr(0.f);
        v3f pe;    // Point on emitter
        v3f ne;    // Surface normal at point
        float pdf; // PDF of choosing point
        sampleEmitterPosition(sampler, em, ne, pe, pdf); // Sample mesh uniformly
        if(scene.bvh->intersect(ray, surfaceInteraction)) {
            v3f emission = getEmission(surfaceInteraction);
            if (emission != v3f(0.f)) {
                return emission;
            }
            v3f wiW = glm::normalize(pe-surfaceInteraction.p);
            surfaceInteraction.wi = surfaceInteraction.frameNs.toLocal(wiW);
            float cosThetaO = glm::dot(ne, -wiW);
            float distance = glm::distance(surfaceInteraction.p, pe);
            Ray shadowRay = Ray(surfaceInteraction.p, wiW);
            double G(0.f);
            v3f brdf = getBSDF(surfaceInteraction)->eval(surfaceInteraction);
            G = glm::max(cosThetaO,0.f) / glm::pow(distance, 2);
            if(m_traceShadows){
                if(scene.bvh->intersect(shadowRay, shadowSurfaceInteracion)){
                    v3f emission = getEmission(shadowSurfaceInteracion);
                    Lr = brdf* emission*G/pdf;
                }
            }else{
                Lr = brdf* G*em.getRadiance()/pdf;
            }
        }



            // TODO(A4): Implement this

        return Lr;
    }

    /// Branch to corresponding method
    v3f render(const Ray& ray, Sampler& sampler) const override {
        switch (m_method) {
            case EPolygonalMethod::ESurfaceArea:
                return PolygonalIntegrator::renderArea(ray, sampler);
                break;
            case EPolygonalMethod::EControlVariates:
                return PolygonalIntegrator::renderControlVariates(ray, sampler);
                break;
            default:
                return PolygonalIntegrator::renderAnalytic(ray, sampler);
                break;
        }
    }

};

TR_NAMESPACE_END