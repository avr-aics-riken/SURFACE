rd /s /q build
cmake -Bbuild -H. -G "Visual Studio 12 Win64" -DSURFACE_BUILD_WITH_OPENMP=On
