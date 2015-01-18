SET(GLES_SOURCES
  gles_accel_builder.cc
  gles_buffer.cc
  gles_c_api.cc
  gles_common.cc
  gles_context.cc
  gles_fragment_op.cc
  gles_glsl_runtime.cc
  gles_graphics_state.cc
  gles_primitive.cc
  gles_program.cc
  gles_raytrace_engine.cc
  gles_render_graph.cc
  gles_resource_manager.cc
  gles_shader.cc
  gles_special.cc
  gles_texture.cc
  gles_vertex.cc

  ../render/accel_line.cc
  ../render/accel_particle.cc
  ../render/accel_tetra.cc
  ../render/accel_triangle.cc
  ../render/camera.cc
  ../render/prefix_tree_util.cc
  ../render/prim_mesh.cc
  ../render/texture.cc
  ../render/toplevel_bvh.cc
  )

IF (SURFACE_WITH_BUILD_MPI)
  ADD_DEFINITIONS(-DLSGL_ENABLE_MPI=1) 
ENDIF()

# Enable UINT64_C
# http://stackoverflow.com/questions/986426/what-do-stdc-limit-macros-and-stdc-constant-macros-mean
ADD_DEFINITIONS(-D__STDC_CONSTANT_MACROS)

INCLUDE_DIRECTORIES(../render ../include)

ADD_LIBRARY(LSGLES ${GLES_SOURCES})
SET_TARGET_PROPERTIES(LSGLES
  PROPERTIES VERSION ${SURFACE_VERSION} SOVERSION ${SURFACE_SOVERSION})

# Install 
INSTALL(TARGETS LSGLES DESTINATION lib)
INSTALL(DIRECTORY ../include/GLES2 DESTINATION include FILES_MATCHING PATTERN "*.h")
INSTALL(DIRECTORY ../include/KHR DESTINATION include FILES_MATCHING PATTERN "*.h")
