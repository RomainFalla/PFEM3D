SET(SRCS
main.cpp
./mesh/Node.hpp ./mesh/Mesh.cpp ./mesh/Mesh.hpp ./mesh/Mesh.inl
./solver/Solver.hpp ./solver/Solver.inl ./solver/Solver.cpp
./quadrature/gausslegendre.hpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})
TARGET_LINK_LIBRARIES(main CGAL::CGAL)
