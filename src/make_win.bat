@echo off
rem This creates the version.h file. 
rem You need git installed (obviously)
rem And to be in the folder where the ".git" directory exists.

FOR /F "delims=" %%i IN ('git describe --tags') DO set GITRESULT=%%i
echo #define ClonalFrameML_GITRevision %GITRESULT% > version.h

rem The linux make.sh file now compiles the code.
rem If you're in VS, remember you need _CRT_SECURE_NO_WARNINGS in the
rem Pre-Processor code.
