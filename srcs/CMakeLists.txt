SET(LIBMESHSRCS
./mesh/Node.hpp ./mesh/Mesh.cpp ./mesh/Mesh.hpp ./mesh/Mesh.inl ./quadrature/gausslegendre.hpp)

SET(LIBSOLVERSRCS
./mesh/Node.hpp ./mesh/Mesh.hpp ./mesh/Mesh.inl ./quadrature/gausslegendre.hpp
./solver/SolverIncompressible.cpp ./solver/SolverIncompressible.hpp ./solver/SolverIncompressible.inl
./solver/SolverCompressible.cpp ./solver/SolverCompressible.hpp ./solver/SolverCompressible.inl)

SET(SRCS
main.cpp
./mesh/Node.hpp ./mesh/Mesh.hpp ./mesh/Mesh.inl ./quadrature/gausslegendre.hpp
./solver/SolverCompressible.hpp ./solver/SolverCompressible.inl
./solver/SolverIncompressible.hpp ./solver/SolverIncompressible.inl)

ADD_LIBRARY(Mesh SHARED ${LIBMESHSRCS})
TARGET_LINK_LIBRARIES(Mesh ${GMSH_LIBRARIES})
TARGET_LINK_LIBRARIES(Mesh CGAL::CGAL)


ADD_LIBRARY(Solver SHARED ${LIBSOLVERSRCS})
TARGET_LINK_LIBRARIES(Solver ${GMSH_LIBRARIES})
TARGET_LINK_LIBRARIES(Solver Mesh)
TARGET_LINK_LIBRARIES(Solver CGAL::CGAL)

ADD_EXECUTABLE(pfem ${SRCS})
TARGET_LINK_LIBRARIES(pfem ${GMSH_LIBRARIES})
TARGET_LINK_LIBRARIES(pfem Mesh)
TARGET_LINK_LIBRARIES(pfem Solver)
TARGET_LINK_LIBRARIES(pfem CGAL::CGAL)
