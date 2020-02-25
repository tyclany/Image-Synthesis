/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/core.h>
#include "core/renderpass.h"
#include "tiny_obj_loader.h"
#include "integrators/path.h"

TR_NAMESPACE_BEGIN

/**
 * Global Illumination baking renderpass.
 */
struct GIPass : RenderPass {
    GLuint shader{0};

    GLuint modelMatUniform{0};
    GLuint viewMatUniform{0};
    GLuint projectionMatUniform{0};

    int m_samplePerVertex;

    std::unique_ptr<PathTracerIntegrator> m_ptIntegrator;

    explicit GIPass(const Scene& scene) : RenderPass(scene) {
        m_ptIntegrator = std::unique_ptr<PathTracerIntegrator>(new PathTracerIntegrator(scene));
        m_ptIntegrator->m_maxDepth = scene.config.integratorSettings.gi.maxDepth;
        m_ptIntegrator->m_rrProb = scene.config.integratorSettings.gi.rrProb;
        m_ptIntegrator->m_rrDepth = scene.config.integratorSettings.gi.rrDepth;
        m_samplePerVertex = scene.config.integratorSettings.gi.samplesByVertex;
    }

    virtual void buildVBO(size_t objectIdx) override {
        GLObject& obj = objects[objectIdx];
        obj.nVerts = scene.getObjectNbVertices(objectIdx);
        for (size_t i = 0; i < (size_t)obj.nVerts; i++) {
            v3f pos = scene.getObjectVertexPosition(objectIdx, i);
            v3f normal = scene.getObjectVertexNormal(objectIdx, i);
            SurfaceInteraction hit = SurfaceInteraction();
            hit.p = pos + normal * Epsilon;
            hit.wo = v3f(1, 2, 3);
            hit.shapeID = objectIdx;
            hit.primID = scene.getPrimitiveID(i);
            hit.frameNg = Frame(normal);
            hit.frameNs = Frame(normal);
            hit.matID = scene.getMaterialID(objectIdx, hit.primID);
            SurfaceInteraction buffer = hit;
            v3f color(0.f);
            Ray ray = Ray(v3f(0.f), v3f(1, 2, 3));
            Sampler sampler = Sampler(260726999);
            for (int j = 0; j < m_samplePerVertex; j++) {
                color += m_ptIntegrator->renderExplicit(ray, sampler,buffer);
            }
            color /= m_samplePerVertex;
            obj.vertices.push_back(pos.x);
            obj.vertices.push_back(pos.y);
            obj.vertices.push_back(pos.z);
            obj.vertices.push_back(color.x);
            obj.vertices.push_back(color.y);
            obj.vertices.push_back(color.z);
            //std::cout << color.x << " " << color.y <<" "<< color.z << endl;
        }
        // TODO(A5): Implement this
        // VBO
        glGenVertexArrays(1, &obj.vao);
        glBindVertexArray(obj.vao);

        glGenBuffers(1, &obj.vbo);
        glBindBuffer(GL_ARRAY_BUFFER, obj.vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     sizeof(GLfloat) * obj.nVerts * N_ATTR_PER_VERT,
                     (GLvoid*) (&obj.vertices[0]),
                     GL_STATIC_DRAW);
    }

    bool init(const Config& config) override {
        RenderPass::init(config);

        // Create shader
        GLuint vs = compileShader("gi.vs", GL_VERTEX_SHADER);
        GLuint fs = compileShader("gi.fs", GL_FRAGMENT_SHADER);
        shader = compileProgram(vs, fs);
        glDeleteShader(vs);
        glDeleteShader(fs);

        // Create uniforms
        modelMatUniform = GLuint(glGetUniformLocation(shader, "model"));
        viewMatUniform = GLuint(glGetUniformLocation(shader, "view"));
        projectionMatUniform = GLuint(glGetUniformLocation(shader, "projection"));

        // Create vertex buffers
        objects.resize(scene.worldData.shapes.size());
        for (size_t i = 0; i < objects.size(); i++) {
            buildVBO(i);
            buildVAO(i);
        }
        return true;
    }

    void cleanUp() override {
        // Delete vertex buffers
        for (size_t i = 0; i < objects.size(); i++) {
            glDeleteBuffers(1, &objects[i].vbo);
            glDeleteVertexArrays(1, &objects[i].vao);
        }

        RenderPass::cleanUp();
    }

    void render() override {
        glBindFramebuffer(GL_FRAMEBUFFER, postprocess_fboScreen);
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        glUseProgram(shader);

        glm::mat4 model, view, projection;
        camera.Update();
        camera.GetMatricies(projection, view, model);

        glUniformMatrix4fv(modelMatUniform, 1, GL_FALSE, &(modelMat[0][0]));
        glUniformMatrix4fv(viewMatUniform, 1, GL_FALSE, &(view[0][0]));
        glUniformMatrix4fv(projectionMatUniform, 1, GL_FALSE, &(projection[0][0]));

        for (auto& object : objects) {
            glBindVertexArray(object.vao);
            glDrawArrays(GL_TRIANGLES, 0, object.nVerts);
            glBindVertexArray(0);
        }
        // TODO(A5): Implement this


        RenderPass::render();
    }

};

TR_NAMESPACE_END
