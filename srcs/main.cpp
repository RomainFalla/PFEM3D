#include <chrono>
#include <iostream>
#include <fstream>

#include "solver/SolverIncompressible.hpp"
#include "solver/SolverCompressible.hpp"

using Clock = std::chrono::high_resolution_clock;
using TimeType = std::chrono::time_point<std::chrono::high_resolution_clock>;

static void displayDT(TimeType startTime, TimeType endTime, std::string text)
{
    auto ellapsedTimeMeasure = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << text << static_cast<double>(ellapsedTimeMeasure.count())/1000.0 << " s" << std::endl;
}

/**
 * \param  argv[1] .json file that contains the parameters.
 * \param  argv[2] .msh file that contains the initial set of nodes.
 */
int main(int argc, char **argv)
{
    //Check that the file format is valid
    if (argc < 3)
    {
        std::cerr   << "Usage: " << argv[0] << " params.json mesh.msh" <<  std::endl;
        return 1;
    }

    try
    {
        std::string paramsFileName {argv[1]};
        std::ifstream paramFile(paramsFileName);
        nlohmann::json j;

        if(!paramFile.is_open())
            throw std::runtime_error("The params file cannot be read!");

        paramFile >> j;
        paramFile.close();

        auto startTime {Clock::now()};

        if(j["ProblemType"] == "Incompressible")
        {
            SolverIncompressible solver(j, std::string(argv[2]));
            solver.solveProblem();
        }
        else if(j["ProblemType"] == "Compressible")
        {
            SolverCompressible solver(j, std::string(argv[2]));
            solver.solveProblem();
        }
        else
            throw std::runtime_error("Unsupported problem type!");

        auto endTime {Clock::now()};
        displayDT(startTime, endTime, "Ellapsed time for problem solving: ");
    }
    catch(const std::bad_alloc& badAllocException)
    {
        std::cerr << "Something went wrong while allocating memory: "
                  << badAllocException.what() << std::endl;

        return -2;
    }
    catch(const std::exception& exception)
    {
        std::cerr << "Something went wrong during program execution: "
                  << exception.what() << std::endl;

        return -1;
    }

    return 0;
}
