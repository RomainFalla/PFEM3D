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
#Please lauch buildCGAL.sh once to build CGAL 5.0 as shared lib

export CGAL_DIR=/c/tools/CGAL-5.0.2
export GMSHSDK=/c/tools/gmsh-4.5.4-Windows64-sdk #put gmsh sdk here
export EIGENSDK=/c/tools/msys64/mingw64/include/eigen3

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lin:${PATH}
export INCLUDE=${GMSHSDK}/include:${EIGENSDK}
export LIB=${GMSHSDK}/lib

cd ../../

rm -rf build
mkdir build
cd build

cmake -G "CodeBlocks - MinGW Makefiles" -DCGAL_HEADER_ONLY=OFF -DUSE_MARCH=1 -DCMAKE_SH=SH-NOTFOUND  ..
cp -r ${GMSHSDK}/lib/gmsh-4.5.dll bin
cp -r ../run/windows/start.bat bin
cp -r /c/tools/msys64/mingw64/bin/libgcc_s_seh-1.dll bin
cp -r /c/tools/msys64/mingw64/bin/libgmp-10.dll bin
cp -r /c/tools/msys64/mingw64/bin/libgmpxx-4.dll bin
cp -r /c/tools/msys64/mingw64/bin/libgomp-1.dll bin
cp -r /c/tools/msys64/mingw64/bin/libstdc++-6.dll bin
cp -r /c/tools/msys64/mingw64/bin/libwinpthread-1.dll bin
cp -r /c/tools/CGAL-5.0.2/bin/libCGAL.dll bin

sed -i 's?C:/building/msys64/mingw64/lib/libmpfr.a?C:/tools/msys64/mingw64/lib/libmpfr.a?g' srcs/CMakeFiles/pfem.dir/build.make
sed -i 's?C:/building/msys64/mingw64/lib/libgmp.dll.a?C:/tools/msys64/mingw64/lib/libgmp.dll.a?g' srcs/CMakeFiles/pfem.dir/build.make
sed -i 's? C:/building/msys64/mingw64/lib/libmpfr.a C:/building/msys64/mingw64/lib/libgmp.dll.a?  C:/tools/msys64/mingw64/lib/libmpfr.a C:/tools/msys64/mingw64/lib/libgmp.dll.a?g' srcs/CMakeFiles/pfem.dir/linklibs.rsp

sed -i 's?C:/building/msys64/mingw64/lib/libmpfr.a?C:/tools/msys64/mingw64/lib/libmpfr.a?g' srcs/CMakeFiles/Mesh.dir/build.make
sed -i 's?C:/building/msys64/mingw64/lib/libgmp.dll.a?C:/tools/msys64/mingw64/lib/libgmp.dll.a?g' srcs/CMakeFiles/Mesh.dir/build.make
sed -i 's? C:/building/msys64/mingw64/lib/libmpfr.a C:/building/msys64/mingw64/lib/libgmp.dll.a?  C:/tools/msys64/mingw64/lib/libmpfr.a C:/tools/msys64/mingw64/lib/libgmp.dll.a?g' srcs/CMakeFiles/Mesh.dir/linklibs.rsp

sed -i 's?C:/building/msys64/mingw64/lib/libmpfr.a?C:/tools/msys64/mingw64/lib/libmpfr.a?g' srcs/CMakeFiles/Solver.dir/build.make
sed -i 's?C:/building/msys64/mingw64/lib/libgmp.dll.a?C:/tools/msys64/mingw64/lib/libgmp.dll.a?g' srcs/CMakeFiles/Solver.dir/build.make
sed -i 's? C:/building/msys64/mingw64/lib/libmpfr.a C:/building/msys64/mingw64/lib/libgmp.dll.a?  C:/tools/msys64/mingw64/lib/libmpfr.a C:/tools/msys64/mingw64/lib/libgmp.dll.a?g' srcs/CMakeFiles/Solver.dir/linklibs.rsp

cd ../
