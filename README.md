[![Build Status](https://travis-ci.org/ImperatorS79/PFEM.svg?branch=master)](https://travis-ci.org/ImperatorS79/PFEM) [![Maintenance](https://img.shields.io/badge/Version-1.0.0-e67e22.svg)](https://github.com/ImperatorS79/PFEM/releases/tag/1.0.0)

# Particular Finite Element Method
This is an implementation of the PFEM method for incompressible and weakly compressible newtonian flow in C++. It uses notably [CGAL](https://www.cgal.org/) for the remeshing, [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) for linear algebra operations and [gmsh](https://www.gmsh.info/) for mesh reading and post-processing.

## Compilation procedure
You should have a C++17 compliant compiler in order to build this project, as well as the different librairies mentionned above. The project also use the [CMake](https://cmake.org/) build system.

### Ubuntu 20.04
You can install the following packages:
```
sudo apt install build-essential libcgal-dev libeigen3-dev liblua5.3-dev swig 4.0 libpython3-dev
```

Then run the file `genCBproject.sh` in the directory `run/linux/` to generate a Code::Blocks project as well as a unix makefile. 

### Windows
The compilation of the project on windows is only supported through MinGW-w64 compiler. It is highly recommanded to use [MSYS2](https://www.msys2.org/) for managing dependencies. You can install the following packages from a MinGW64 shell:

```
pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-toolchain mingw-w64-x86_64-boost mingw-w64-x86_64-mpfr mingw-w64-x86_64-gmp mingw-w64-x86_64-nlohmann-json mingw-w64-x86_64-eigen3 mingw-w64-x86_64-lua mingw-w64-x86_64-swig
```

Then run the file `genCBprojectMSYS2.sh` in the directory `run/windows/` (inside the MinGW64 shell) to generate a Code::Blocks project as well as a mingw makefile. 

### MacOS
The compilation of the project on macos is only supported through homebrew's clang (as it needs openmp). You can install the following packages from homebrew:

```
brew update && brew upgrade cmake && brew install wget gnu-tar llvm libomp unzip swig cgal lua python@3.7 nlohmann-json eigen;
```

Then run the file `genCBproject.sh` in the directory `run/macos/` to generate a Code::Blocks project as well as a unix makefile.
