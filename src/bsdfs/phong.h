/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Modified Phong reflectance model
 */
struct PhongBSDF : BSDF {

    std::unique_ptr<Texture < v3f>> specularReflectance;
    std::unique_ptr<Texture < v3f>> diffuseReflectance;
    std::unique_ptr<Texture < float>> exponent;
    float specularSamplingWeight;
    float scale;

    PhongBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.specular_texname.empty())
            specularReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.specular)));
        else
            specularReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.specular_texname));

        if (mat.diffuse_texname.empty())
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        exponent = std::unique_ptr<Texture<float>>(new ConstantTexture1f(mat.shininess));

        //get scale value to ensure energy conservation
        v3f maxValue = specularReflectance->getMax() + diffuseReflectance->getMax();
        float actualMax = max(max(maxValue.x, maxValue.y), maxValue.z);
        scale = actualMax > 1.0f ? 0.99f * (1.0f / actualMax) : 1.0f;

        float dAvg = getLuminance(diffuseReflectance->getAverage() * scale);
        float sAvg = getLuminance(specularReflectance->getAverage() * scale);
        specularSamplingWeight = sAvg / (dAvg + sAvg);

        components.push_back(EGlossyReflection);
        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (unsigned int component : components)
            combinedType |= component;
    }

    inline float getExponent(const SurfaceInteraction& i) const override {
        return exponent->eval(worldData, i);
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    v3f eval(const SurfaceInteraction& i) const override {
        v3f val(0.f);
        if(Frame::cosTheta(i.wi)>0&&Frame::cosTheta(i.wo)>0) {
            v3f diffuseReflectanceBuffer = (PhongBSDF::scale)*diffuseReflectance->eval(worldData, i);
            v3f specularReflectanceBuffer = (PhongBSDF::scale)*specularReflectance->eval(worldData,i);
            float exponentBuffer = exponent->eval(worldData,i);
            v3f rl = PhongBSDF::reflect(i.wi);
            rl=normalize(rl);
            // normalize to get the correct cosine value
            float cosineValue=glm::dot(rl,i.wo);
            // cosine is the angle between the reflected wave and the incident wave
            float foreShortening = Frame::cosTheta(i.wi);
            //handout instruction
            val=foreShortening*((diffuseReflectanceBuffer*INV_PI)+specularReflectanceBuffer*(exponentBuffer+2)*INV_TWOPI*glm::max(0.f,(glm::pow(cosineValue,exponentBuffer))));
        }
        // TODO(A2): Implement this

        return val;
    }

    float pdf(const SurfaceInteraction& i) const override {
        float pdf = 0.f;
        float pdfPhong = 0.f;
        v3f wiDiffuse = i.wi;
        float pdfDiffuse=Warp::squareToCosineHemispherePdf(wiDiffuse);
        v3f wr = reflect(i.wo);
        Frame reflectFrame = Frame(wr);
        wiDiffuse = reflectFrame.toLocal(wiDiffuse);
        if(wiDiffuse.z<0){
            pdfPhong = 0;
        }else{
            pdfPhong=Warp::squareToPhongLobePdf(wiDiffuse,exponent->eval(worldData,i));
        }
        // TODO(A3): Implement this
        pdf = pdfPhong*specularSamplingWeight+pdfDiffuse*(1-specularSamplingWeight);
        return pdf;
    }

    v3f sample(SurfaceInteraction& i, Sampler& sampler, float* pdfMix) const override {
        v3f val(0.f);
        if(sampler.next()<specularSamplingWeight){
            v3f wiPhong = Warp::squareToPhongLobe(sampler.next2D(), exponent->eval(worldData, i));
            v3f wr = reflect(i.wo);
            Frame reflectFrame = Frame(wr);
            i.wi = reflectFrame.toWorld(wiPhong);
            *pdfMix = pdf(i);
            val = eval(i);
        }
        else{
            v3f wiDiffuse = Warp::squareToCosineHemisphere(sampler.next2D());
            i.wi = wiDiffuse;
            *pdfMix = pdf(i);
            val = eval(i);
        }
        // TODO(A3): Implement this

        return val;
    }

    std::string toString() const override { return "Phong"; }
};

TR_NAMESPACE_END