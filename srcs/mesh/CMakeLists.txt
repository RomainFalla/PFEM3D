find_path(EIGEN_INCLUDE_DIRS NAMES "Eigen")
if(NOT EIGEN_INCLUDE_DIRS)
    message(FATAL_ERROR "Eigen not found!")
else()
	message(STATUS "Found Eigen: " ${EIGEN_INCLUDE_DIRS})
endif()

find_package(CGAL REQUIRED)
if(CGAL_FOUND)
    message(STATUS "Found CGAL: " ${CGAL_DIR})
endif()

set(LIBMESHSRCS
./Node.hpp ./Node.inl
./Element.hpp ./Element.inl ./Element.cpp
./Facet.hpp ./Facet.inl ./Facet.cpp
./Mesh.hpp ./Mesh.inl ./Mesh.cpp ./Mesh2D.cpp ./Mesh3D.cpp
./mesh_defines.h ../sharedLib_defines.h)

add_library(pfemMesh SHARED ${LIBMESHSRCS})
target_include_directories(pfemMesh SYSTEM
                           PRIVATE ${EIGEN_INCLUDE_DIRS}
                           PRIVATE ${GMSH_INCLUDE_DIRS})
target_link_libraries(pfemMesh PRIVATE OpenMP::OpenMP_CXX CGAL::CGAL ${GMSH_LIBRARIES})
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    target_compile_options(pfemMesh PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES CLANG)
    target_compile_options(pfemMesh PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wnull-dereference -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    target_compile_options(pfemMesh PRIVATE /W4 /WX /wd4251)
endif()
