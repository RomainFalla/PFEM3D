#include <chrono>
#include <exception>
#include <iostream>
#include <string>
#include <fstream>

#include <nlohmann/json.hpp>

#include "solver/SolverIncompressible.hpp"
#include "solver/SolverCompressible.hpp"

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
        std::string paramsFileName = std::string(argv[1]);
        std::ifstream paramFile(paramsFileName);
        nlohmann::json j;

        if(!paramFile.is_open())
            throw std::runtime_error("The params file cannot be read!");

        paramFile >> j;
        paramFile.close();

        auto startTime = std::chrono::high_resolution_clock::now();

        if(j["ProblemType"] == "Incompressible")
        {
            SolverIncompressible solver(j, std::string(argv[2]), std::string(argv[3]));
            solver.solveProblem();
        }
        else if(j["ProblemType"] == "Compressible")
        {
            SolverCompressible solver(j, std::string(argv[2]), std::string(argv[3]));
            solver.solveProblem();
        }
        else
            throw std::runtime_error("Unsupported problem type!");

        auto endTime = std::chrono::high_resolution_clock::now();
        auto ellapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
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
