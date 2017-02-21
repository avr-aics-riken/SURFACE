SET(SRCS
  main.cc
  )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
ADD_DEFINITIONS(-D__STDC_CONSTANT_MACROS)

ADD_EXECUTABLE(variable_line_size_render ${SRCS})
ADD_SANITIZERS(variable_line_size_render)

IF (WIN32)
TARGET_LINK_LIBRARIES(variable_line_size_render LSGLES)
ELSE()
TARGET_LINK_LIBRARIES(variable_line_size_render LSGLES dl)
ENDIF()

# --
INSTALL(TARGETS variable_line_size_render DESTINATION bin)
INSTALL(FILES input.frag DESTINATION bin) # fixme
# --
