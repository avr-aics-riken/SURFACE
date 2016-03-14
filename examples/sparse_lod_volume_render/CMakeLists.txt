SET(SRCS
  main.cc
  )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
ADD_DEFINITIONS(-D__STDC_CONSTANT_MACROS)

ADD_EXECUTABLE(sparse_lod_volume_render ${SRCS})

IF (WIN32)
TARGET_LINK_LIBRARIES(sparse_lod_volume_render LSGLES)
ELSE()
TARGET_LINK_LIBRARIES(sparse_lod_volume_render LSGLES dl)
ENDIF()

# --
INSTALL(TARGETS sparse_lod_volume_render DESTINATION bin)
# --
