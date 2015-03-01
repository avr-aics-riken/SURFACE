#include "../../gles/gles_c_api.h"
#include <GLES2/gl2ext.h>

#include "../common/SimpleTGA.h"

int kWindowWidth = 512;
int kWindowHeight = 512;

namespace {

bool SaveColorBufferRGBA(const char *savename) {
  void *tgabuffer;
  unsigned char *imgBuf = new unsigned char[kWindowWidth * kWindowHeight * 4];
  glReadPixels(0, 0, kWindowWidth, kWindowHeight, GL_RGBA, GL_UNSIGNED_BYTE,
               imgBuf);
  int tgasize =
      SimpleTGASaverRGBA(&tgabuffer, kWindowWidth, kWindowHeight, imgBuf);
  delete[] imgBuf;
  if (!tgasize) {
    printf("Failed save.\n");
    return false;
  }

  FILE *fp = fopen(savename, "wb");
  fwrite(tgabuffer, 1, tgasize, fp);
  fclose(fp);
  free(tgabuffer);

  return true;
}

bool LoadShader(GLuint &prog, GLuint &fragShader,
                       const char *fragShaderFilename) {
  GLint val = 0;

  // free old shader/program
  if (prog != 0)
    glDeleteProgram(prog);
  if (fragShader != 0)
    glDeleteShader(fragShader);

  static GLchar srcbuf[16384];
  FILE *fp = fopen(fragShaderFilename, "rb");
  if (!fp) {
    return false;
  }

  fseek(fp, 0, SEEK_END);
  size_t len = ftell(fp);
  rewind(fp);
  len = fread(srcbuf, 1, len, fp);
  srcbuf[len] = 0;
  fclose(fp);

  const GLchar *src = srcbuf;

  fragShader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragShader, 1, &src, NULL);
  glCompileShader(fragShader);
  glGetShaderiv(fragShader, GL_COMPILE_STATUS, &val);
  assert(val == GL_TRUE && "failed to compile shader");

  prog = glCreateProgram();
  glAttachShader(prog, fragShader);
  glLinkProgram(prog);
  glGetProgramiv(prog, GL_LINK_STATUS, &val);
  assert(val == GL_TRUE && "failed to link shader");

  return true;
}

}  // namespace

int main(int argc, char *argv[]) {
  const char *frag_shader_file_name = NULL;

  if (argc < 2) {
    frag_shader_file_name = "default.frag";
  } else {
    frag_shader_file_name = argv[1];
  }

  GLuint prog = 0, frag_shader = 0;
  if (!LoadShader(prog, frag_shader , frag_shader_file_name)) {
    fprintf(stderr, "failed to load shader: %s\n", frag_shader_file_name);
    return -1;
  }

  assert(glIsProgram(prog) == GL_TRUE);

  glUseProgram(prog);

  // Create frame buffer
  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // Create color render buffer
  GLuint color_renderbuffer;
  glGenRenderbuffers(1, &color_renderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, color_renderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA, kWindowWidth, kWindowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                            GL_RENDERBUFFER, color_renderbuffer);

  // Create depth render buffer
  GLuint depth_renderbuffer;
  glGenRenderbuffers(1, &depth_renderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, depth_renderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT,
                        kWindowWidth, kWindowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                            GL_RENDERBUFFER, depth_renderbuffer);

  // Clear background color
  glClearColor(0, 0, 0, 1);
  glClear(GL_COLOR_BUFFER_BIT);

  // Set viewport
  glViewport(0, 0, kWindowWidth, kWindowHeight);

  // Set camera
  GLfloat eye[3] = {0.0, 0.0, 2.0};
  GLfloat lookat[3] = {0.0, 0.0, 0.0};
  GLfloat up[3] = {0.0, 1.0, 0.0};

  lsglSetCamera(eye, lookat, up, 45.0f);

  // Set parameters to the shader
  glUniform1f(glGetUniformLocation(prog, "time"), 100.0);
  glUniform2f(glGetUniformLocation(prog, "resolution"),
              static_cast<float>(kWindowWidth),
              static_cast<float>(kWindowHeight));

  // Create vertex buffer

  GLfloat vertex_data[] = {
    -1.0f, -1.0f, 0.0f,
     1.0f, -1.0f, 0.0f,
    -1.0f,  1.0f, 0.0f,

    -1.0f,  1.0f, 0.0f,
     1.0f, -1.0f, 0.0f,
     1.0f,  1.0f, 0.0f
  };

  GLuint vertex_buffer;
  glGenBuffers(1, &vertex_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
  lsglBufferDataPointer(GL_ARRAY_BUFFER, sizeof(vertex_data),
                        vertex_data, GL_STATIC_DRAW);

  // Set shader attribute
  GLint attr_ptr = glGetAttribLocation(prog, "position");
  glVertexAttribPointer(attr_ptr, 3, GL_FLOAT, GL_FALSE,
                        sizeof(float) * 3, static_cast<void *>(0));
  glEnableVertexAttribArray(attr_ptr);
  assert(glGetError() == GL_NO_ERROR);

  // Draw
  glDrawArrays(GL_TRIANGLES, 0, sizeof(vertex_data) / sizeof(vertex_data[0]));

  glFinish();

  // Save the image
  const char *save_image_file_name = "colorbuf.tga";
  if (!SaveColorBufferRGBA(save_image_file_name)) {
    fprintf(stderr, "failed to save the image to file: %s\n",
            save_image_file_name);
    return 1;
  }

  glDeleteRenderbuffers(1, &color_renderbuffer);
  glDeleteRenderbuffers(1, &depth_renderbuffer);
  glDeleteFramebuffers(1, &framebuffer);

  return 0;
}
