[![Build Status](https://travis-ci.org/ImperatorS79/PFEM.svg?branch=master)](https://travis-ci.org/ImperatorS79/PFEM) [![Maintenance](https://img.shields.io/badge/Version-1.0.0-e67e22.svg)](https://github.com/ImperatorS79/PFEM/releases/tag/1.0.0)

# Particular Finite Element Method
This is an implementation of the PFEM method for incompressible and weakly compressible newtonian flow in C++. It uses notably [CGAL](https://www.cgal.org/) for the remeshing, [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) for linear algebra operations and [gmsh](https://www.gmsh.info/) for mesh reading and post-processing.

## Compilation procedure
You should have a C++14 compliant compiler in order to build this project, as well as the different librairies mentionned above. The project also use the [CMake](https://cmake.org/) build system.

### Ubuntu 19.10
You can install the following packages:
```
sudo apt install build-essential libcgal-dev libeigen3-dev liblua5.3-dev
```

Do not forget to download gmsh (it will be downloaded automatically if you use the `run/linux/genCBproject.sh` file). 

### Windows
It is highly recommended to build the project using MinGW-w64 compiler from [MSYS2](https://www.msys2.org/). You can install the following packages from a MinGW64 shell:

```
pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-toolchain mingw-w64-x86_64-eigen3 mingw-w64-x86_64-nlohmann-json mingw-w64-x86_64-lua
```

You should also download [CGAL](https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-5.0.2/CGAL-5.0.2-Setup.exe) and install it. Then run the `buildCGAL.sh` file in a MinGW64 shell.
You should also download [sol3/sol.hpp](https://github.com/ThePhD/sol2/releases/download/v3.0.3/sol.hpp) and [sol3/forward.hpp](https://github.com/ThePhD/sol2/releases/download/v3.0.3/forward.hpp).

The cmake command should be append by `-DCMAKE_SH=SH-NOTFOUND` if you use cmake from the MinGW64 shell (the easiest way to have everything in `PATH`)

Do not forget to download gmsh!

### MacOS
Nothing should prevent the project to be built on macOS and it is nearly like linux (it is highly recommende to install a compiler supporting OpenMP).
