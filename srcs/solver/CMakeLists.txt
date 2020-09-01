FIND_PATH(EIGEN_INCLUDE_DIRS NAMES "Eigen")
if(NOT EIGEN_INCLUDE_DIRS)
    MESSAGE(FATAL_ERROR "Eigen not found!")
ELSE()
	MESSAGE(STATUS "Found Eigen: " ${EIGEN_INCLUDE_DIRS})
ENDIF()

FIND_PACKAGE(Lua 5.1 REQUIRED)

FIND_PATH(SOL3_INCLUDE_DIRS NAMES "sol3/sol.hpp")
if(NOT SOL3_INCLUDE_DIRS)
    MESSAGE(FATAL_ERROR "sol3 not found!")
ELSE()
	MESSAGE(STATUS "Found sol3: " ${SOL3_INCLUDE_DIRS})
ENDIF()

SET(LIBSOLVERSRCS
./extractors/Extractor.cpp ./extractors/PointExtractor.cpp
./extractors/GMSHExtractor.cpp ./extractors/MinMaxExtractor.cpp ./extractors/MassExtractor.cpp
./extractors/Extractor.hpp ./extractors/PointExtractor.hpp ./extractors/MassExtractor.hpp
./extractors/GMSHExtractor.hpp ./extractors/MinMaxExtractor.hpp
./Solver.cpp
./Solver.hpp ./Solver.inl
./SolverIncompressible.cpp
./SolverIncompressible.hpp
./SolverCompressible.cpp
./SolverCompressible.hpp ./SolverCompressible.inl)

ADD_LIBRARY(pfemSolver SHARED ${LIBSOLVERSRCS})
GENERATE_EXPORT_HEADER(pfemSolver
                       EXPORT_MACRO_NAME SOLVER_API
					   EXPORT_FILE_NAME Solver_export.h)
TARGET_INCLUDE_DIRECTORIES(pfemSolver SYSTEM
                           PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
                           PRIVATE ${GMSH_INCLUDE_DIRS}
                           PUBLIC ${EIGEN_INCLUDE_DIRS}
                           PUBLIC ${LUA_INCLUDE_DIR}
                           PUBLIC ${SOL3_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(pfemSolver
                      PRIVATE ${GMSH_LIBRARIES}
                      PUBLIC ${LUA_LIBRARIES} pfemMesh)
IF(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    TARGET_COMPILE_OPTIONS(pfemSolver PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wshadow)
ELSEIF(CMAKE_CXX_COMPILER_ID MATCHES CLANG)
    TARGET_COMPILE_OPTIONS(pfemSolver PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wnull-dereference -Wshadow)
ELSEIF(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    TARGET_COMPILE_OPTIONS(pfemSolver PRIVATE /W4 /WX /wd4251)
ENDIF()

FIND_PACKAGE(Python 3 REQUIRED COMPONENTS Interpreter Development)
MESSAGE(STATUS "Found python: " ${Python_INCLUDE_DIRS})
MESSAGE(STATUS "Found python: " ${Python_LIBRARY_DIRS})
MESSAGE(STATUS "Found python: " ${Python_RUNTIME_LIBRARY_DIRS})
cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)
FIND_PACKAGE(SWIG 4.0 REQUIRED)
INCLUDE(${SWIG_USE_FILE})

IF(PYTHONLIBS_VERSION_STRING MATCHES "^2.*$" )
   set(CMAKE_SWIG_FLAGS -classic)
ELSE()
   set(CMAKE_SWIG_FLAGS -py3)
ENDIF()

SET(CMAKE_SWIG_OUTDIR ${PROJECT_BINARY_DIR}/bin/ CACHE PATH "")

SET_PROPERTY(SOURCE ./Solver.i PROPERTY CPLUSPLUS ON)
SWIG_ADD_LIBRARY(pfemSolverw
				 TYPE SHARED
				 LANGUAGE python
				 SOURCES ./Solver.i)
TARGET_INCLUDE_DIRECTORIES(pfemSolverw PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
TARGET_INCLUDE_DIRECTORIES(pfemSolverw PUBLIC ${Python_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
SWIG_LINK_LIBRARIES(pfemSolverw
                   pfemSolver ${Python_LIBRARIES})
