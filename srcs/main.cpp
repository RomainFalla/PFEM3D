#include <chrono>
#include <iostream>
#include <string>
#if defined(_OPENMP)
    #include <cstdlib>
    #include <omp.h>
#endif

#include "mesh/Node.hpp"
#include "mesh/Mesh.hpp"
#include "params/Params.hpp"
#include "solver/Solver.hpp"

/**
 * @param  argv[1] .msh file that contains the mesh.
 * @param  argv[2] .dat file that contains the parameters.
 * @param  argv[3] name of the .msh file that will contain the results.
 */
int main(int argc, char **argv)
{
    // check that the file format is valid
    if (argc < 2)
    {
        std::cerr   << "Usage: " << argv[0] << " param.json"
                    <<  std::endl;
        return 1;
    }

    Params params;
    if(!params.loadFromFile(std::string(argv[1])))
    {
        std::cerr   << "Something went wrong when trying to read the file" << argv[1]
                    <<  std::endl;
        return 1;
    }

	Mesh mesh(params);
	mesh.loadFromFile("dummyName");

    Solver solver(params, mesh, INCOMPRESSIBLE_PSPG);
    solver.solveProblem();
	
    return 0;
}
