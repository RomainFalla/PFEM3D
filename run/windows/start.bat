::Copy this file in the build/Release/bin or build/Debug/bin directory

@ECHO OFF
SET OMP_NUM_THREADS=2
pfem.exe params/params.json geometry/square.msh results.msh
PAUSE