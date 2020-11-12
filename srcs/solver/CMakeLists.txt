find_path(EIGEN_INCLUDE_DIRS NAMES "Eigen")
if(NOT EIGEN_INCLUDE_DIRS)
    message(FATAL_ERROR "Eigen not found!")
else()
	message(STATUS "Found Eigen: " ${EIGEN_INCLUDE_DIRS})
endif()

find_package(Lua 5.1 REQUIRED)

find_path(SOL3_INCLUDE_DIRS NAMES "sol3/sol.hpp")
if(NOT SOL3_INCLUDE_DIRS)
    message(FATAL_ERROR "sol3 not found!")
else()
	message(STATUS "Found sol3: " ${SOL3_INCLUDE_DIRS})
endif()

set(LIBSOLVERSRCS
./extractors/Extractor.cpp ./extractors/PointExtractor.cpp
./extractors/GMSHExtractor.cpp ./extractors/MinMaxExtractor.cpp ./extractors/MassExtractor.cpp
./extractors/Extractor.hpp ./extractors/PointExtractor.hpp ./extractors/MassExtractor.hpp
./extractors/GMSHExtractor.hpp ./extractors/MinMaxExtractor.hpp
./Solver.cpp ./Solver.hpp ./Solver.inl
./SolverIncompressible.cpp ./SolverIncompressible.hpp
./SolverCompressible.cpp ./SolverCompressible.hpp ./SolverCompressible.inl
./SolverConduction.cpp ./SolverConduction.hpp
./SolverBoussinesq.cpp ./SolverBoussinesq.hpp
./SolverBoussinesqWC.cpp ./SolverBoussinesqWC.hpp
./solver_defines.h ../sharedLib_defines.h)

add_library(pfemSolver SHARED ${LIBSOLVERSRCS})
target_include_directories(pfemSolver SYSTEM
                           PRIVATE ${GMSH_INCLUDE_DIRS}
                           PUBLIC ${EIGEN_INCLUDE_DIRS}
                           PUBLIC ${LUA_INCLUDE_DIR}
                           PUBLIC ${SOL3_INCLUDE_DIRS})
target_link_libraries(pfemSolver
                      PRIVATE ${GMSH_LIBRARIES}
                      PUBLIC OpenMP::OpenMP_CXX
                      PUBLIC ${LUA_LIBRARIES} pfemMesh)
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    target_compile_options(pfemSolver PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES CLANG)
    target_compile_options(pfemSolver PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wnull-dereference -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    target_compile_options(pfemSolver PRIVATE /W4 /WX /wd4251)
endif()

find_package(Python 3 REQUIRED COMPONENTS Interpreter Development)
message(STATUS "Found python: " ${Python_INCLUDE_DIRS})
message(STATUS "Found python: " ${Python_LIBRARY_DIRS})
message(STATUS "Found python: " ${Python_RUNTIME_LIBRARY_DIRS})
cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)
find_package(SWIG 4.0 REQUIRED)
INCLUDE(${SWIG_USE_FILE})

if(PYTHONLIBS_VERSION_STRING MATCHES "^2.*$" )
   set(CMAKE_SWIG_FLAGS -classic)
else()
   set(CMAKE_SWIG_FLAGS -py3)
endif()

set(CMAKE_SWIG_OUTDIR ${PROJECT_BINARY_DIR}/bin/ CACHE PATH "")

set_property(SOURCE ./Solver.i PROPERTY CPLUSPLUS ON)
swig_add_library(pfemSolverw
				 TYPE SHARED
				 LANGUAGE python
				 SOURCES ./Solver.i)
target_include_directories(pfemSolverw PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(pfemSolverw PUBLIC ${Python_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
swig_link_libraries(pfemSolverw
                   pfemSolver ${Python_LIBRARIES})
