#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>

#include "simulation/physics/Problems.hpp"
#include "simulation/utility/Utility.hpp"
#include "simulation/utility/SolTable.hpp"

/**
 * \param  argv[1] .lua file that contains the parameters.
 */
int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr   << "Usage: " << argv[0] << " params.lua" <<  std::endl;
        return 1;
    }

    std::unique_ptr<Problem> pProblem;

    try
    {
        std::ifstream luaFile(argv[1]);
        if(!luaFile.is_open())
            throw std::runtime_error("Cannot open lua file " + std::string(argv[1]) + "!");

        luaFile.close();

        //hack
        sol::state state;
        state.script_file(argv[1]);

        SolTable table = SolTable("Problem", state);

        std::string problemType = table.checkAndGet<std::string>("id");

        if(problemType == "IncompNewtonNoT" || problemType == "Boussinesq" || problemType == "Conduction")
            pProblem = std::make_unique<ProbIncompNewton>(argv[1]);
        else if(problemType == "WCompNewtonNoT" || problemType == "BoussinesqWC")
            pProblem = std::make_unique<ProbWCompNewton>(argv[1]);
        else
            throw std::runtime_error("Unknown problem type " + problemType + "!");

        pProblem->simulate();
    }
    catch(const std::exception& e)
    {
        std::cerr << "Something went wrong while running the program: " << e.what() << std::endl;
        return -1;
    }
    catch(...)
    {
        std::cerr << "An unknown exception has occurred." << std::endl;
        return -1;
    }

    return 0;
}
