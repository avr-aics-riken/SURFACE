GLES2_INCS=-I../../include
GLES2_LIBS=-L../../gles -lLSGLES
CXX=g++
CXXFLAGS=-O2 -g -D__STDC_CONSTANT_MACROS
EXTRA_LIBS=-ldl
# for windows
#GLES2_LIBS=-L../../gles ../../gles/LSGLES.lib -lwinmm
#CXX=x86_64-w64-mingw32-g++.exe

all:
	$(CXX) $(GLES2_INCS) $(CXXFLAGS) main.cc tinymt64.cpp $(GLES2_LIBS) $(EXTRA_LIBS)
