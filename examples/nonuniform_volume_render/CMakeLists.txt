SET(SRCS
  main.cc
  )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
ADD_DEFINITIONS(-D__STDC_CONSTANT_MACROS)

ADD_EXECUTABLE(nonuniform_volume_render ${SRCS})
ADD_SANITIZERS(nonuniform_volume_render)

IF (WIN32)
TARGET_LINK_LIBRARIES(nonuniform_volume_render LSGLES)
ELSE()
TARGET_LINK_LIBRARIES(nonuniform_volume_render LSGLES dl)
ENDIF()

# --
INSTALL(TARGETS nonuniform_volume_render DESTINATION bin)
# --
