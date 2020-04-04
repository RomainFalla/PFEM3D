::Copy this file in the build/bin directory

@ECHO OFF
SET OMP_NUM_THREADS=2
pfem.exe ../../params/xD/params.json ../../geometry/xD/mesh.msh results.msh
PAUSE