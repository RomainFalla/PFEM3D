#!/bin/sh
#Example of usage of PFEM
#Please take care of the lib location !

export GMSHSDK=../../dependencies/gmsh-4.6.0-MacOSX-sdk/

export OMP_NUM_THREADS=2

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:"${PATH}"

./pfem ../../examples/2D/pipe/paramsIncomp.json ../../examples/2D/pipe/geometry.msh
