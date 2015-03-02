#ifdef LSGL_ENABLE_MPI
#include <mpi.h>
#endif

#include "../../gles/gles_c_api.h"
#include <GLES2/gl2ext.h>

#include "../common/SimpleTGA.h"

#include "easywsclient.hpp"
#include "jpge.h"

#ifdef _WIN32
#pragma comment( lib, "ws2_32" )
#include <WinSock2.h>
#endif

int kWindowWidth = 512;
int kWindowHeight = 512;
const char* kHost = "localhost";
int kPort = 8080;

namespace {

using easywsclient::WebSocket;
static WebSocket::pointer ws = NULL;

void handle_message(const std::string & message)
{
    printf(">>> %s\n", message.c_str());
    if (message == "ack") { ws->close(); }
}

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

static unsigned char*
EncodeAsJPEG(int& out_size, unsigned char *rgba, int w, int h)
{
    size_t sz = w * h * 4;  // alloc enough size.
    if (sz <  1024) sz = 1024; // at lest 1K is required.
    out_size = sz;


    //   num_channels must be 1 (Y), 3 (RGB), or 4 (RGBA), image pitch must be width*num_channels.
    //   bool compress_image_to_jpeg_file(const char *pFilename, int width, int height, int num_channels, 
    //                                    const uint8 *pImage_data, const params &comp_params = params());
    jpge::params comp_params = jpge::params();
    comp_params.m_quality = 100;
    void *buf = malloc(sz);
    bool ret = jpge::compress_image_to_jpeg_file_in_memory(buf, out_size, w, h, 4, rgba, comp_params);

    return reinterpret_cast<unsigned char*>(buf); // callee must call free() 
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

  int rank = 0;
#ifdef LSGL_ENABLE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

#ifdef _WIN32
    INT rc;
    WSADATA wsaData;

    rc = WSAStartup(MAKEWORD(2, 2), &wsaData);
    if (rc) {
        printf("WSAStartup Failed.\n");
        return 1;
    }
#endif

  const char *frag_shader_file_name = NULL;

  if (argc < 2) {
    frag_shader_file_name = "default.frag";
  } else {
    frag_shader_file_name = argv[1];
  }

  if (argc > 2) {
    kHost = argv[2];
  }

  if (argc > 3) {
    kPort = atoi(argv[3]);
  }

  
  char buf[8192];
  sprintf(buf, "ws://%s:%d/", kHost, kPort);
  std::string addr(buf);

  if (rank == 0) {
    ws = WebSocket::from_url(addr);
    if (!ws) {
      // fail to connect server. save image as jpg as file.
    } else {
      //ws->send("goodbye");
      //ws->send("hello");
      //while (ws->getReadyState() != WebSocket::CLOSED) {
      //  ws->poll();
      //  ws->dispatch(handle_message);
      //}
      //delete ws;
    }
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
  glDrawArrays(GL_TRIANGLES, 0, sizeof(vertex_data) / sizeof(vertex_data[0]) / 3);

  glFinish();

  if (rank == 0 && ws) {

    int dataLen = 0;
    unsigned char *imgBuf = new unsigned char[kWindowWidth * kWindowHeight * 4];
    glReadPixels(0, 0, kWindowWidth, kWindowHeight, GL_RGBA, GL_UNSIGNED_BYTE, imgBuf);
    unsigned char* data = EncodeAsJPEG(dataLen, imgBuf, kWindowWidth, kWindowHeight);

    std::string s = std::string(reinterpret_cast<const char*>(data), reinterpret_cast<const char*>(data) + dataLen);

    printf("dataLen = %d\n", dataLen);
    ws->sendBinary(s);
    //while (ws->getReadyState() != WebSocket::CLOSED) {
    //  ws->poll();
    //  ws->dispatch(handle_message);
    //}

    while (ws->getReadyState() != WebSocket::CLOSED) {
      ws->poll();
      ws->dispatch(handle_message);
    }

    delete [] imgBuf;
    free(data);

  }

  // Save the image
  char namebuf[1024];
#ifdef LSGL_ENABLE_MPI
  sprintf(namebuf, "colorbuf_%04d.tga", rank);
#else
  sprintf(namebuf, "colorbuf.tga");
#endif
  if (!SaveColorBufferRGBA(namebuf)) {
    fprintf(stderr, "failed to save the image to file: %s\n",
            namebuf);
    return 1;
  }

  glDeleteRenderbuffers(1, &color_renderbuffer);
  glDeleteRenderbuffers(1, &depth_renderbuffer);
  glDeleteFramebuffers(1, &framebuffer);

#ifdef LSG_ENABLE_MPI
  MPI_Finalize();
#endif

#ifdef _WIN32
  WSACleanup();
#endif

  return 0;
}
