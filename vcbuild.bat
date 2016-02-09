chcp 437
call "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\vcvarsall.bat" x86_amd64
rem TODO: build all projects.
msbuild build/SURFACE.sln /p:Configuration=Release /t:LSGLES:Rebuild
