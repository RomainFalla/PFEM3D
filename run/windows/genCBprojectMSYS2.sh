#!/bin/sh
#Launch this in MinGW64 shell from MSYS2
#You should have installed the following package from MSYS2:
#mingw-w64-x86_64-cmake
#mingw-w64-x86_64-toolchain
#mingw-w64-x86_64-eigen3
#mingw-w64-x86_64-nlohmann-json
#https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-5.0.2/CGAL-5.0.2-Setup.exe (CGAL .cmake from MSYS2 bugged)
#Please use buildCGAL.sh If you want to use cgal as a shared lib (recommanded for build time)

export HERE=$PWD

cd /c/tools

export HERE_TOOLS=$PWD
if [ ! -d "cgal-releases-CGAL-5.0.3" ]; then
    wget https://github.com/CGAL/cgal/archive/releases/CGAL-5.0.3.tar.gz
    tar -xf CGAL-5.0.3.tar.gz
    rm -rf CGAL-5.0.3.tar.gz
    cd $HERE_TOOLS
fi

IsShared="static"
if [ ! -z "$1" ]; then
	IsShared="$1"
fi

if [ "$IsShared" != "static" ] && [ "$IsShared" != "shared" ]; then
    echo "Unknown lib type"
	echo "$IsShared"
	exit 1
fi

buildType="Release"
if [ ! -z "$2" ]; then
	buildType="$2"
fi

if [ "$buildType" != "Release" ] && [ "$buildType" != "Debug" ]; then
    echo "Unknown build type"
	echo "$buildTYpe"
	exit 1
fi

if [ "$IsShared" = "shared" ] && [ ! -d  "/c/tools/CGAL-5.0.3/" ]; then
	cd /c/tools/cgal-releases-CGAL-5.0.3
	rm -rf build
	mkdir build
	cd build
	cmake -G "MinGW Makefiles" -DCGAL_HEADER_ONLY=OFF -DCMAKE_INSTALL_PREFIX="/c/tools/CGAL-5.0.3/" -DCMAKE_BUILD_TYPE="${buildType}" .. 
	mingw32-make
	mingw32-make install
    cd $HERE_TOOLS
fi

if [ ! -d "sol3" ]; then
  mkdir sol3
  cd sol3
  wget https://github.com/ThePhD/sol2/releases/download/v3.0.3/sol.hpp
  wget https://github.com/ThePhD/sol2/releases/download/v3.0.3/forward.hpp
  cd ../
fi

if [ ! -d "gmsh-4.6.0-Windows64-sdk" ]; then
  wget http://gmsh.info/bin/Windows/gmsh-4.6.0-Windows64-sdk.zip
  unzip gmsh-4.6.0-Windows64-sdk.zip 
  rm -rf gmsh-4.6.0-Windows64-sdk.zip
fi

if [ "$IsShared" = "shared" ]; then 
	export CGAL_DIR=/c/tools/CGAL-5.0.3
else
	export CGAL_DIR=/c/tools/cgal-releases-CGAL-5.0.3
fi
export GMSHSDK=/c/tools/gmsh-4.6.0-Windows64-sdk #put gmsh sdk here
export EIGENSDK=/c/tools/msys64/mingw64/include/eigen3
export SOL3SDK=/c/tools/

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lin:${PATH}
export INCLUDE=${GMSHSDK}/include:${EIGENSDK}:${SOL3SDK}
export LIB=${GMSHSDK}/lib

cd $HERE

cd ../../

rm -rf build
mkdir build
cd build

if [ "$IsShared" = "shared" ]; then 
	cmake -G "CodeBlocks - MinGW Makefiles" -DCGAL_HEADER_ONLY=OFF -DUSE_MARCH=0 -DCMAKE_SH=SH-NOTFOUND  ..
else
	cmake -G "CodeBlocks - MinGW Makefiles" -DCGAL_HEADER_ONLY=ON -DUSE_MARCH=0 -DCMAKE_SH=SH-NOTFOUND  ..
fi

cp -r ${GMSHSDK}/lib/gmsh-4.6.dll bin
cp -r ../run/windows/start.bat bin
cp -r /c/tools/msys64/mingw64/bin/libgcc_s_seh-1.dll bin
cp -r /c/tools/msys64/mingw64/bin/libgmp-10.dll bin
cp -r /c/tools/msys64/mingw64/bin/libgmpxx-4.dll bin
cp -r /c/tools/msys64/mingw64/bin/libgomp-1.dll bin
cp -r /c/tools/msys64/mingw64/bin/libstdc++-6.dll bin
cp -r /c/tools/msys64/mingw64/bin/libwinpthread-1.dll bin
if [ "$IsShared" = "shared" ]; then 
	cp -r /c/tools/CGAL-5.0.3/bin/libCGAL.dll bin
fi

cp -r /c/tools/msys64/mingw64/bin/lua53.dll bin

cd ../
