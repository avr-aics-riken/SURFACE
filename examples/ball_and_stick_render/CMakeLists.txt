SET(SRCS
  main.cc
  tiny_pdb.cc
  )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
ADD_DEFINITIONS(-D__STDC_CONSTANT_MACROS)

ADD_EXECUTABLE(ball_and_stick_render ${SRCS})

IF (WIN32)
TARGET_LINK_LIBRARIES(ball_and_stick_render LSGLES)
ELSE()
TARGET_LINK_LIBRARIES(ball_and_stick_render LSGLES dl)
ENDIF()

# --
INSTALL(TARGETS ball_and_stick_render DESTINATION bin)
# --
