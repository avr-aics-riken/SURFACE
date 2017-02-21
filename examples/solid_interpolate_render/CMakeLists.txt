SET(SRCS
  main.cc
  tinymt64.cpp
  )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
ADD_DEFINITIONS(-D__STDC_CONSTANT_MACROS)

ADD_EXECUTABLE(solid_interpolate_render ${SRCS})
ADD_SANITIZERS(solid_interpolate_render)

IF (WIN32)
TARGET_LINK_LIBRARIES(solid_interpolate_render LSGLES)
ELSE()
TARGET_LINK_LIBRARIES(solid_interpolate_render LSGLES dl)
ENDIF()

# --
INSTALL(TARGETS solid_interpolate_render DESTINATION bin)
# --
