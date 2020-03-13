SET(LIBMESHSRCS
./mesh/Node.hpp ./mesh/Mesh.cpp ./mesh/Mesh.hpp ./mesh/Mesh.inl ./quadrature/gausslegendre.hpp)

SET(LIBSOLVERSRCS
./mesh/Node.hpp ./mesh/Mesh.hpp ./mesh/Mesh.inl ./quadrature/gausslegendre.hpp
./solver/Solver.hpp ./solver/Solver.inl ./solver/Solver.cpp)

SET(SRCS
main.cpp
./mesh/Node.hpp ./mesh/Mesh.hpp ./mesh/Mesh.inl ./quadrature/gausslegendre.hpp
./solver/Solver.hpp ./solver/Solver.inl
)

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
