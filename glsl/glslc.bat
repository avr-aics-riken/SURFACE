@echo off
set GLSLC_ERRORLEVEL=

setlocal
set path=%~dp0;%~dp0..;%path%
set scriptname=%~dp0glslc
python "%scriptname%" %*
endlocal & set GLSLC_ERRORLEVEL=%ERRORLEVEL%

if NOT "%COMSPEC%" == "%SystemRoot%\system32\cmd.exe" goto returncode
if errorlevel 9009 echo you do not have python in your PATH
goto endglslc

:returncode
exit /B %GLSLC_ERRORLEVEL%

:endglslc
call :returncode %GLSLC_ERRORLEVEL%

