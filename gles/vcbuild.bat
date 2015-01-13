rem -- config --
rem set target=vs2012
rem set buildtype=release

set target=vs2013
set buildtype=release
set myplatform=x64

..\\tools\\windows\\premake5.exe %target% %buildtype%
