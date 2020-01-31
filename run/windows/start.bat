::Copy this file in the build/Release/bin or build/Debug/bin directory

@ECHO OFF
SET OMP_NUM_THREADS=1
SET OMP_CANCELLATION=true
main.exe "geometry/young/young.msh" "params/young.json"
PAUSE