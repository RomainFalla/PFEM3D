# Particular Finite Element Method
This is an implementation of the PFEM method for incompressible and weakly compressible newtonian flow in C++. It uses notably [CGAL](https://www.cgal.org/) for the remeshing, [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) for linear algebra operations and [gmsh](https://www.gmsh.info/) for mesh reading and post-processing.

## Compilation procedure
You should have a C++14 compliant compiler in order to build this project, as well as the different librairies mentionned above. The project also use the [CMake](https://cmake.org/) build system.

### Ubuntu 19.10
You can install the following packages:
```
sudo apt install build-essential libcgal-dev libeigen3-dev
```

Do not forget to download gmsh (it will be downloaded automatically if you use the `run/linux/genCBproject.sh` file).

### Windows
It is highly recommended to build the project using MinGW-w64 compiler from [MSYS2](https://www.msys2.org/). You can install the following packages from a MinGW64 shell:

```
pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-toolchain mingw-w64-x86_64-eigen3 mingw-w64-x86_64-nlohmann-json
wget http://repo.msys2.org/mingw/x86_64/mingw-w64-x86_64-cgal-4.11-1-any.pkg.tar.xz
pacman -U path/to/mingw-w64-x86_64-cgal-4.11-1-any.pkg.tar.xz
```

The cmake command should be append by `-DCMAKE_SH=SH-NOTFOUND` if you use cmake from the MinGW64 shell (the easiest way to have everything in `PATH`)

An older version of CGAL is downloaded due to a bug in the latest version proposed by MSYS2. However another bug affect the chosen version. The files `srcs/CMakeFiles/pfem.dir/build.make`, `srcs/CMakeFiles/pfem.dir/linklibs.rsp`, `srcs/CMakeFiles/Mesh.dir/build.make`, `srcs/CMakeFiles/Mesh.dir/linklibs.rsp`, `srcs/CMakeFiles/Solver.dir/build.make`and `srcs/CMakeFiles/Solver.dir/linklibs.rsp` will contain reference to a `C:/building/msys64` directory. Replace that directory by the installation directory of MSYS2 (everything is handled if you use the `run/windows/genCBproject.sh` file. 

Do not forget to download gmsh (it will be downloaded automatically if you use the `run/windows/genCBproject.sh` file).

### MacOS
Nothing should prevent the project to be built on macOS and it is nearly like linux (it is highly recommende to install a compiler supporting OpenMP).
