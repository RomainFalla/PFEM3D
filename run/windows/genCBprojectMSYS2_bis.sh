cd ../../
if [ ! -d "dependencies" ]; then
    mkdir dependencies
fi

cd dependencies/

if [ ! -d "gmsh-4.7.0-Windows64-sdk" ]; then
  wget http://gmsh.info/bin/Windows/gmsh-4.7.0-Windows64-sdk.zip
  unzip gmsh-4.7.0-Windows64-sdk.zip 
  rm -rf gmsh-4.7.0-Windows64-sdk.zip
fi

export CGAL_DIR=/c/msys64/mingw64/lib/cmake/CGAL/
export GMSHSDK=${PWD}/gmsh-4.7.0-Windows64-sdk #put gmsh sdk here
export EIGENSDK=/c/msys64/mingw64/include/eigen3
export SOL3SDK=/c/mingw64/include/

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lin:${PATH}
export INCLUDE=${GMSHSDK}/include:${EIGENSDK}:${SOL3SDK}
export LIB=${GMSHSDK}/lib

cd ../

rm -rf build
mkdir build
cd build

cmake -G "CodeBlocks - MinGW Makefiles" -DCGAL_HEADER_ONLY=ON -DUSE_MARCH=0 -DCMAKE_SH=SH-NOTFOUND  ..

cp -r ${GMSHSDK}/lib/gmsh-4.7.dll bin
cp -r ../run/windows/start.bat bin
cp -r /c/msys64/mingw64/bin/libgcc_s_seh-1.dll bin
cp -r /c/msys64/mingw64/bin/libgmp-10.dll bin
cp -r /c/msys64/mingw64/bin/libgmpxx-4.dll bin
cp -r /c/msys64/mingw64/bin/libgomp-1.dll bin
cp -r /c/msys64/mingw64/bin/libstdc++-6.dll bin
cp -r /c/msys64/mingw64/bin/libwinpthread-1.dll bin
cp -r /c/msys64/mingw64/bin/lua53.dll bin

cd ../
