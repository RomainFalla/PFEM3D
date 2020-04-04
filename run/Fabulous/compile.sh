#!/bin/bash
#Do not forget to install libgfortran3  on your system!

module load gcc/6.3.0 boost/1.78.0 gmp/6.0.2 mpfr/4.0.2 cmake/current CGAL/4.14.3 eigen/3.3.7

cd ../../

if [ ! -d "gmsh-4.5.4-Linux64-sdk" ]; then
  wget http://gmsh.info/bin/Linux/gmsh-4.5.4-Linux64-sdk.tgz
  tar -xf gmsh-4.5.4-Linux64-sdk.tgz 
  rm -rf gmsh-4.5.4-Linux64-sdk.tgz 
fi

if [ ! -d "nlohmann" ]; then
  wget https://github.com/nlohmann/json/releases/download/v3.6.1/include.zip
  unzip include.zip
  rm -rf include.zip
  mkdir nlohmann/
  mv include/nlohmann nlohmann/nlohmann
  rm -rf include
fi

export GMSHSDK=${PWD}/gmsh-4.5.4-Linux64-sdk/
export JSONSDK=${PWD}/nlohmann/
export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:"${PATH}"
export INCLUDE=${GMSHSDK}/include:"${INCLUDE}"
export INCLUDE=${JSONSDK}:"${INCLUDE}"
export LIB=${GMSHSDK}/lib:"${LIB}"
export PYTHONPATH=${GMSHSDK}/lib:"${PYTHONPATH}"
export DYLD_LIBRARY_PATH=${GMSHSDK}/lib:"${DYLD_LIBRARY_PATH}"

rm -rf build
mkdir build
cd build

cmake ../ -DCMAKE_BUILD_TYPE=Release -DUSE_MARCH=1  -G "Unix Makefiles"

make

cd ../