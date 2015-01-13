#include "gtest/gtest.h"
#include "gles_context.h"
#include "gles_resource_manager.h"

using namespace lsgl;

TEST(VertexAttribTest, BufferData) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      

  const float cube[] = {  -1.0, -1.0, -1.0,
                          -1.0, -1.0,  1.0,
                          -1.0,  1.0, -1.0,
                          -1.0,  1.0,  1.0,
                           1.0, -1.0, -1.0,
                           1.0, -1.0,  1.0,
                           1.0,  1.0, -1.0,
                           1.0,  1.0,  1.0 };

  const GLuint cube_indices[] = {
      0, 2, 1, 1, 2, 3,
      4, 5, 6, 5, 7, 6,
      0, 1, 4, 1, 5, 4,
      2, 7, 3, 2, 6, 7,
      1, 3, 5, 3, 7, 5,
      0, 4, 2, 2, 4, 6 };

  // Create Vertex Buffers.
  GLuint cubevtx, cubeidx;
  ctx.glGenBuffers(1, &cubevtx);
  ctx.glBindBuffer(GL_ARRAY_BUFFER, cubevtx);
  ctx.glBufferData(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  ctx.glGenBuffers(1, &cubeidx);
  ctx.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeidx);
  ctx.glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cube_indices), cube_indices, GL_STATIC_DRAW);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  glBindBuffer(GL_ARRAY_BUFFER, cubevtx);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeidx);
  GLuint pos = 0;
  glVertexAttribPointer(pos, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, (void*)0);
  glEnableVertexAttribArray(pos);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
}

