#!/bin/sh
#Launch this in MinGW64 shell from MSYS2
#You should have installed the following package from MSYS2:
#mingw-w64-x86_64-cmake
#mingw-w64-x86_64-toolchain
#mingw-w64-x86_64-eigen3
#mingw-w64-x86_64-nlohmann-json
#https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-5.0.2/CGAL-5.0.2-Setup.exe (CGAL .cmake from MSYS2 bugged)
#After creating the build file you should run "grep -r C:/building" in the build directory and
#replace every "C:/building" by by the directory where msys64 directory is

rm -rf /c/tools/CGAL-5.0.2
buildType="Release"
if [ ! -z "$1" ]; then
	buildType="$1"
fi
cd /c/dev/CGAL-5.0.2
rm -rf build
mkdir build
cd build
cmake -G "MinGW Makefiles" -DCGAL_HEADER_ONLY=OFF -DCMAKE_INSTALL_PREFIX="/c/tools/CGAL-5.0.2/" -DCMAKE_BUILD_TYPE="${buildType}" .. 
mingw32-make
mingw32-make install