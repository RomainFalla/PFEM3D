include(GenerateExportHeader)

SET(LIBMESHSRCS
./mesh/Mesh.cpp ./mesh/Mesh2D.cpp ./mesh/Mesh3D.cpp
./mesh/Mesh.hpp ./mesh/Mesh.inl ./mesh/Node.hpp)

SET(LIBSOLVERSRCS
./solver/extractors/Extractor.cpp ./solver/extractors/PointExtractor.cpp
./solver/extractors/GMSHExtractor.cpp ./solver/extractors/MinMaxExtractor.cpp ./solver/extractors/MassExtractor.cpp
./solver/extractors/Extractor.hpp ./solver/extractors/PointExtractor.hpp ./solver/extractors/MassExtractor.hpp
./solver/extractors/GMSHExtractor.hpp ./solver/extractors/MinMaxExtractor.hpp
./solver/Solver.cpp
./solver/Solver.hpp ./solver/Solver.inl
./solver/SolverIncompressible.cpp
./solver/SolverIncompressible.hpp
./solver/SolverCompressible.cpp
./solver/SolverCompressible.hpp ./solver/SolverCompressible.inl)


SET(SRCS
main.cpp)

ADD_LIBRARY(Mesh SHARED ${LIBMESHSRCS})
GENERATE_EXPORT_HEADER(Mesh
                       EXPORT_MACRO_NAME MESH_API
					   EXPORT_FILE_NAME Mesh_export.h)
TARGET_INCLUDE_DIRECTORIES(Mesh PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
TARGET_LINK_LIBRARIES(Mesh ${GMSH_LIBRARIES})
TARGET_LINK_LIBRARIES(Mesh CGAL::CGAL)

ADD_LIBRARY(Solver SHARED ${LIBSOLVERSRCS})
GENERATE_EXPORT_HEADER(Solver
                       EXPORT_MACRO_NAME SOLVER_API
					   EXPORT_FILE_NAME Solver_export.h)
TARGET_INCLUDE_DIRECTORIES(Solver PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
TARGET_LINK_LIBRARIES(Solver ${GMSH_LIBRARIES})
TARGET_LINK_LIBRARIES(Solver Mesh)

ADD_EXECUTABLE(pfem ${SRCS})
TARGET_INCLUDE_DIRECTORIES(pfem PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
TARGET_LINK_LIBRARIES(pfem ${GMSH_LIBRARIES})
TARGET_LINK_LIBRARIES(pfem Mesh)
TARGET_LINK_LIBRARIES(pfem Solver)
