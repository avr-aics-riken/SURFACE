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

  ../render/render_accel_line.cc
  ../render/render_accel_particle.cc
  ../render/render_accel_solid.cc
  ../render/render_accel_tetra.cc
  ../render/render_accel_triangle.cc
  ../render/render_accel_volume.cc
  ../render/render_camera.cc
  ../render/render_prefix_tree_util.cc
  ../render/render_prim_triangle.cc
  ../render/render_texture.cc
  ../render/render_bvh_tree.cc
  ../render/tinymt64.cpp
  )

SET(LSGLES_SUFFIX "")
SET(DEBUG_PREFIX "")

IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  SET(DEBUG_PREFIX "d")
ENDIF()

IF (SURFACE_BUILD_WITH_MPI)
  IF (SURFACE_BUILD_SCREEEN_PARALLEL)
    ADD_DEFINITIONS(-DLSGL_ENABLE_SCREEN_PARALLEL=1) 
    SET(LSGLES_SUFFIX "${DEBUG_PREFIX}_mpi_screenparallel")
  ELSE()
    SET(LSGLES_SUFFIX "${DEBUG_PREFIX}_mpi")
  ENDIF()
ELSE()
  SET(LSGLES_SUFFIX "${DEBUG_PREFIX}")
ENDIF()

STRING(COMPARE NOTEQUAL "${LSGLES_SUFFIX}" "" hassuffix)
IF (hassuffix)
  IF (BUILD_SHARED_LIBS)
    SET(LSGLES_SUFFIX "${LSGLES_SUFFIX}${CMAKE_SHARED_LIBRARY_SUFFIX}")
  ELSE()
    SET(LSGLES_SUFFIX "${LSGLES_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  ENDIF()
ENDIF()

# Enable UINT64_C
# http://stackoverflow.com/questions/986426/what-do-stdc-limit-macros-and-stdc-constant-macros-mean
ADD_DEFINITIONS(-D__STDC_CONSTANT_MACROS)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/render ${PROJECT_SOURCE_DIR}/include)

ADD_LIBRARY(LSGLES ${GLES_SOURCES})
ADD_SANITIZERS(LSGLES)

STRING(COMPARE NOTEQUAL "${LSGLES_SUFFIX}" "" suffixresult)
IF (suffixresult)
  SET_TARGET_PROPERTIES(LSGLES PROPERTIES SUFFIX ${LSGLES_SUFFIX})
ENDIF()
SET_TARGET_PROPERTIES(LSGLES
  PROPERTIES VERSION ${SURFACE_VERSION} SOVERSION ${SURFACE_SOVERSION})

# Install 
INSTALL(TARGETS LSGLES DESTINATION lib)
INSTALL(DIRECTORY ../include/GLES2 DESTINATION include FILES_MATCHING PATTERN "*.h")
INSTALL(DIRECTORY ../include/KHR DESTINATION include FILES_MATCHING PATTERN "*.h")
