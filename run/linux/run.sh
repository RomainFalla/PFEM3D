#!/bin/sh
#Example of usage of PFEM
#Please take care of the lib location !


export GMSHSDK=../../../gmsh-4.5.2-Linux64-sdk/

export OMP_NUM_THREADS=2

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:"${PATH}"

./main params.json square.msh results.msh
