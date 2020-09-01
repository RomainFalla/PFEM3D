#!/bin/sh
#Do not forget to install libgfortran3  on your system!
#You will also need need lua 5.3 (i.e. liblua5.3-dev), CGAL 4.14 (i.e. libcgal-dev), swig 4.0 (swig4.0) and python 3 (i.e. libpython3-dev) (at least)
#This project needs a C++17 compliant compiler

cd ../../
if [ ! -d "dependencies" ]; then
    mkdir dependencies
fi

cd dependencies/

if [ ! -d "sol3" ]; then
  mkdir sol3
  cd sol3
  wget https://github.com/ThePhD/sol2/releases/download/v3.0.3/sol.hpp
  wget https://github.com/ThePhD/sol2/releases/download/v3.0.3/forward.hpp
  cd ../
fi
if [ ! -d "gmsh-4.6.0-Linux64-sdk" ]; then
  wget http://gmsh.info/bin/Linux/gmsh-4.6.0-Linux64-sdk.tgz
  tar -xf gmsh-4.6.0-Linux64-sdk.tgz 
  rm -rf gmsh-4.6.0-Linux64-sdk.tgz 
fi

if [ ! -d "eigen-3.3.7" ]; then
  wget https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
  tar -xf eigen-3.3.7.tar.gz
  rm -rf eigen-3.3.7.tar.gz
fi

if [ ! -d "nlohmann" ]; then
  wget https://github.com/nlohmann/json/releases/download/v3.9.1/include.zip
  unzip include.zip
  rm -rf include.zip
  mkdir nlohmann/
  mv include/nlohmann nlohmann/nlohmann
  rm -rf include
fi

export GMSHSDK=${PWD}/gmsh-4.6.0-Linux64-sdk/
export EIGENSDK=${PWD}/eigen-3.3.7/
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

cd ../

rm -rf build
mkdir build
cd build

cmake ../ -DCMAKE_BUILD_TYPE=Debug  -G "CodeBlocks - Unix Makefiles"

cd ../
