@rem change to build directory
cd /d %~dp0
@rem call "C:\Program Files\Microsoft Visual Studio 10.0\VC\bin\vcvars32.bat"

nmake -y -f nmake.mak tidyall
