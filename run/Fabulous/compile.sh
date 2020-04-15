#!/bin/bash
#Do not forget to install libgfortran3  on your system!

module load gcc/9.3.0 boost/1.78.0 gmp/6.0.2 mpfr/4.0.2 cmake/current CGAL/4.14.3 eigen/3.3.7

cd ../../

if [ ! -d "lua" ]; then
  wget https://www.lua.org/ftp/lua-5.3.5.tar.gz
  tar -xf lua-5.3.5.tar.gz 
  rm -rf lua-5.3.5.tar.gz
  cd lua-5.3.5
  sed -i 's?INSTALL_TOP= /usr/local?INSTALL_TOP= ../../lua?g' Makefile
  sed -i 's?CFLAGS= -O2 -Wall -Wextra -DLUA_COMPAT_5_2 $(SYSCFLAGS) $(MYCFLAGS)?CFLAGS= -O2 -Wall -Wextra -fPIC -DLUA_COMPAT_5_2 $(SYSCFLAGS) $(MYCFLAGS)?g' src/Makefile
  
  make linux
  make install
  cd ../
  rm -rf lua-5.3.5 
fi

if [ ! -d "sol3" ]; then
  mkdir sol3
  cd sol3
  wget https://github.com/ThePhD/sol2/releases/download/v3.0.3/sol.hpp
  wget https://github.com/ThePhD/sol2/releases/download/v3.0.3/forward.hpp
  cd ../
fi

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
export SOLSDK=${PWD}
export LUASDK=${PWD}/lua/
export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:${LUASDK}/bin:${LUASDK}/lib"${PATH}"
export INCLUDE=${GMSHSDK}/include:${LUASDK}/include:"${INCLUDE}"
export INCLUDE=${JSONSDK}:"${INCLUDE}"
export INCLUDE=${SOLSDK}:"${INCLUDE}"
export LIB=${GMSHSDK}/lib:${LUASDK}/lib:"${LIB}"
export PYTHONPATH=${GMSHSDK}/lib:"${PYTHONPATH}"
export DYLD_LIBRARY_PATH=${GMSHSDK}/lib:"${DYLD_LIBRARY_PATH}"

rm -rf build
mkdir build
cd build

cmake ../ -DCMAKE_BUILD_TYPE=Release -DUSE_MARCH=1  -G "Unix Makefiles"

make

cd ../