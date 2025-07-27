#include "overlay.h"
#include <iostream>
#include <cstring> // For memcpy

RenderObject* createRenderObject(GLenum mode, const float obj_color[3], float obj_size) {
    RenderObject* object = new RenderObject();
    if (!object) {
        std::cerr << "Failed to allocate memory for RenderObject." << std::endl;
        return nullptr;
    }

    object->drawing_mode = mode;
    object->size = obj_size;
    memcpy(object->color, obj_color, 3 * sizeof(float));

    glGenVertexArrays(1, &object->vao);
    glGenBuffers(1, &object->vbo);

    glBindVertexArray(object->vao);
    glBindBuffer(GL_ARRAY_BUFFER, object->vbo);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);

    return object;
}

void destroyRenderObject(RenderObject* object) {
    if (!object) return;
    glDeleteVertexArrays(1, &object->vao);
    glDeleteBuffers(1, &object->vbo);
    delete object;
}

void updateRenderObject(RenderObject* object) {
    if (!object) return;
    // It's okay to update an object with zero vertices (clearing it)
    if (object->vertices.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, object->vbo);
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        return;
    }
    glBindBuffer(GL_ARRAY_BUFFER, object->vbo);
    glBufferData(GL_ARRAY_BUFFER, object->vertices.size() * sizeof(Vertex), object->vertices.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void drawRenderObject(const RenderObject* object, unsigned int shader_program_id) {
    if (!object || object->vertices.empty()) return;

    // Set the per-object properties before drawing
    // Note: The shader program is assumed to be in use already by the caller.

    // 1. Set the color uniform
    GLint color_loc = glGetUniformLocation(shader_program_id, "objectColor");
    if (color_loc != -1) {
        glUniform3fv(color_loc, 1, object->color);
    }

    // 2. Set the size (for points or lines)
    if (object->drawing_mode == GL_POINTS) {
        glPointSize(object->size);
    } else if (object->drawing_mode == GL_LINES || object->drawing_mode == GL_LINE_STRIP || object->drawing_mode == GL_LINE_LOOP) {
        glLineWidth(object->size);
    }

    // 3. Bind the VAO and draw the object
    glBindVertexArray(object->vao);
    glDrawArrays(object->drawing_mode, 0, (GLsizei)object->vertices.size());
    glBindVertexArray(0);
}