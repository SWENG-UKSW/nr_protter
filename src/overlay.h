#ifndef RENDER_OBJECTS_H
#define RENDER_OBJECTS_H

#include <vector>
#include <glad/glad.h>

// Represents a single 2D vertex
struct Vertex {
    float x, y;
};

// Expanded structure to hold all necessary data for a renderable object.
struct RenderObject {
    unsigned int vao = 0;
    unsigned int vbo = 0;
    GLenum drawing_mode;
    std::vector<Vertex> vertices;
    float color[3]; // RGB color
    float size;     // For point size or line width
};

/**
* @brief Creates a new RenderObject with a specific color and size.
*
* @param mode The OpenGL drawing mode (e.g., GL_POINTS, GL_LINES).
* @param obj_color A 3-element float array for the object's RGB color.
* @param obj_size The size (for points) or width (for lines).
* @return A pointer to the newly created RenderObject.
*/
RenderObject* createRenderObject(GLenum mode, const float obj_color[3], float obj_size);

/**
* @brief Frees the memory and OpenGL resources used by a RenderObject.
*
* @param object The RenderObject to destroy.
*/
void destroyRenderObject(RenderObject* object);

/**
* @brief Updates the object's vertex buffer on the GPU.
*
* @param object The RenderObject to update.
*/
void updateRenderObject(RenderObject* object);

/**
* @brief Draws the RenderObject, setting its unique color and size uniforms.
*
* @param object The RenderObject to draw.
* @param shader_program_id The ID of the shader program to use for drawing.
*/
void drawRenderObject(const RenderObject* object, unsigned int shader_program_id);

#endif // RENDER_OBJECTS_H