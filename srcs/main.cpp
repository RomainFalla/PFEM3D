#include <chrono>
#include <exception>
#include <iostream>
#include <string>
#include <fstream>
#if defined(_OPENMP)
    #include <cstdlib>
    #include <omp.h>
#endif

#include "mesh/Node.hpp"
#include "mesh/Mesh.hpp"
#include "params/Params.hpp"
#include "solver/Solver.hpp"

/**
 * @param  argv[1] .json file that contains the parameters.
 * @param  argv[2] .msh file that contains the initial set of nodes.
 * @param  argv[3] name of the .msh file that will contain the results.
 */
int main(int argc, char **argv)
{
    // check that the file format is valid
    if (argc < 4)
    {
        std::cerr   << "Usage: " << argv[0] << " params.json mesh.msh results.msh"
                    <<  std::endl;
        return 1;
    }

    try
    {
        Params params;
        params.loadFromFile(std::string(argv[1]));

        auto startTime = std::chrono::high_resolution_clock::now();
        Mesh mesh(params);
        mesh.loadFromFile(std::string(argv[2]));
        auto endTime = std::chrono::high_resolution_clock::now();
        auto ellapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        std::cout << "Ellapsed time for mesh loading: "
                  << static_cast<double>(ellapsedTime.count())/1000.0
                  << " s" << std::endl;

        startTime = std::chrono::high_resolution_clock::now();
        Solver solver(params, mesh, std::string(argv[3]));
        solver.solveProblem();
        endTime = std::chrono::high_resolution_clock::now();
        ellapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        std::cout << "Ellapsed time for problem solving: "
                  << static_cast<double>(ellapsedTime.count())/1000.0
                  << " s" << std::endl;
    }
    catch(std::exception const & exception)
    {
        std::cerr << "Something went wrong during program execution: "
                  << exception.what() << std::endl;

        return -1;
    }

    return 0;
}
