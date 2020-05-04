#!/bin/sh
#Do not forget to install libgfortran3  on your system!
#You will also need need lua 5.3 and CGAL 4.14 (at least)
#This project needs a C++17 compliant compiler

cd ../../
if [ ! -d "sol3" ]; then
  mkdir sol3
  cd sol3
  wget https://github.com/ThePhD/sol2/releases/download/v3.0.3/sol.hpp
  wget https://github.com/ThePhD/sol2/releases/download/v3.0.3/forward.hpp
  cd ../
fi
if [ ! -d "gmsh-4.5.6-Linux64-sdk" ]; then
  wget http://gmsh.info/bin/Linux/gmsh-4.5.6-Linux64-sdk.tgz
  tar -xf gmsh-4.5.6-Linux64-sdk.tgz 
  rm -rf gmsh-4.5.6-Linux64-sdk.tgz 
fi

if [ ! -d "eigen-eigen-323c052e1731" ]; then
  wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
  tar -xf 3.3.7.tar.gz
  rm -rf 3.3.7.tar.gz
fi

if [ ! -d "nlohmann" ]; then
  wget https://github.com/nlohmann/json/releases/download/v3.6.1/include.zip
  unzip include.zip
  rm -rf include.zip
  mkdir nlohmann/
  mv include/nlohmann nlohmann/nlohmann
  rm -rf include
fi

export GMSHSDK=${PWD}/gmsh-4.5.6-Linux64-sdk/
export EIGENSDK=${PWD}/eigen-eigen-323c052e1731/
export JSONSDK=${PWD}/nlohmann/
export SOLSDK=${PWD}/

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:"${PATH}"
export INCLUDE=${GMSHSDK}/include:"${INCLUDE}"
export INCLUDE=${EIGENSDK}:"${INCLUDE}"
export INCLUDE=${JSONSDK}:"${INCLUDE}"
export INCLUDE=${SOLSDK}:"${INCLUDE}"
export LIB=${GMSHSDK}/lib:"${LIB}"
export PYTHONPATH=${GMSHSDK}/lib:"${PYTHONPATH}"
export DYLD_LIBRARY_PATH=${GMSHSDK}/lib:"${DYLD_LIBRARY_PATH}"

rm -rf build
mkdir build
cd build

cmake ../ -DCMAKE_BUILD_TYPE=Debug  -G "CodeBlocks - Unix Makefiles"

cd ../
