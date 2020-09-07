#include <chrono>
#include <iostream>
#include <fstream>

#include <nlohmann/json.hpp>

#include "solver/SolverIncompressible.hpp"
#include "solver/SolverCompressible.hpp"


using Clock = std::chrono::high_resolution_clock;
using TimeType = std::chrono::time_point<std::chrono::high_resolution_clock>;

static void displayDT(TimeType startTime, TimeType endTime, std::string text)
{
    auto ellapsedTimeMeasure = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << text << static_cast<double>(ellapsedTimeMeasure.count())/1000.0 << " s" << std::endl;
}

static void addExctractors(Solver& solver, nlohmann::json j)
{
    auto extractors = j["Solver"]["Extractors"];
    for(auto& extractor : extractors)
    {
        if(extractor["type"].get<std::string>() == "Point")
        {
            solver.addPointExtractor(extractor["outputFile"].get<std::string>(),
                                     extractor["timeBetweenWriting"].get<double>(),
                                     extractor["stateToWrite"].get<unsigned short>(),
                                     extractor["points"].get<std::vector<std::vector<double>>>());
        }
        else if(extractor["type"].get<std::string>() == "GMSH")
        {
            solver.addGMSHExtractor(extractor["outputFile"].get<std::string>(),
                                    extractor["timeBetweenWriting"].get<double>(),
                                    extractor["whatToWrite"].get<std::vector<std::string>>(),
                                    extractor["writeAs"].get<std::string>());
        }
        else if(extractor["type"].get<std::string>() == "MinMax")
        {
            solver.addMinMaxExtractor(extractor["outputFile"].get<std::string>(),
                                      extractor["timeBetweenWriting"].get<double>(),
                                      extractor["coordinate"].get<unsigned short>(),
                                      extractor["minMax"].get<std::string>());
        }
        else if(extractor["type"].get<std::string>() == "Mass")
        {
            solver.addMassExtractor(extractor["outputFile"].get<std::string>(),
                                    extractor["timeBetweenWriting"].get<double>());
        }
        else
        {
            std::string errotText = std::string("unknown extractor type ") + extractor["Type"].get<std::string>();
            throw std::runtime_error(errotText);
        }
    }
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

        MeshCreateInfo meshInfos = {};
        meshInfos.hchar = j["Remeshing"]["hchar"].get<double>();
        meshInfos.alpha = j["Remeshing"]["alpha"].get<double>();
        meshInfos.omega = j["Remeshing"]["omega"].get<double>();
        meshInfos.gamma = j["Remeshing"]["gamma"].get<double>();
        meshInfos.boundingBox = j["Remeshing"]["boundingBox"].get<std::vector<double>>();
        meshInfos.mshFile = std::string(argv[2]);

        SolverCreateInfo solverInfos = {};
        solverInfos.gravity = j["Solver"]["gravity"].get<double>();
        solverInfos.strongPAtFS = j["Solver"]["strongPAtFS"].get<bool>();
        solverInfos.adaptDT = j["Solver"]["Time"]["adaptDT"].get<bool>();
        solverInfos.initialDT = j["Solver"]["Time"]["initialDT"].get<double>();
        solverInfos.endTime = j["Solver"]["Time"]["endTime"].get<double>();
        solverInfos.maxDT = j["Solver"]["Time"]["maxDT"].get<double>();
        solverInfos.IBCfile = j["Solver"]["IBCs"].get<std::string>();
        solverInfos.meshInfos = std::move(meshInfos);

        bool verboseOutput = j["verboseOutput"].get<bool>();

        auto startTime {Clock::now()};

        if(j["ProblemType"] == "Incompressible")
        {
            SolverIncompCreateInfo solverIncompInfos;
            solverIncompInfos.rho = j["Solver"]["Fluid"]["rho"].get<double>();
            solverIncompInfos.mu = j["Solver"]["Fluid"]["mu"].get<double>();
            solverIncompInfos.picardMaxIter = j["Solver"]["Picard"]["maxIter"].get<double>();
            solverIncompInfos.picardRelTol = j["Solver"]["Picard"]["relTol"].get<double>();
            solverIncompInfos.coeffDTdecrease = j["Solver"]["Time"]["coeffDTdecrease"].get<double>();
            solverIncompInfos.coeffDTincrease = j["Solver"]["Time"]["coeffDTincrease"].get<double>();
            solverIncompInfos.solverInfos = std::move(solverInfos);

            SolverIncompressible solver(std::move(solverIncompInfos));
            addExctractors(solver, j);
            j.clear();
            solver.solveProblem(verboseOutput);
        }
        else if(j["ProblemType"] == "Compressible")
        {
            SolverCompCreateInfo solverCompInfos;
            solverCompInfos.strongContinuity = j["Solver"]["strongContinuity"].get<bool>();
            solverCompInfos.securityCoeff = j["Solver"]["Time"]["securityCoeff"].get<double>();
            solverCompInfos.rho0 = j["Solver"]["Fluid"]["rho0"].get<double>();
            solverCompInfos.mu = j["Solver"]["Fluid"]["mu"].get<double>();
            solverCompInfos.K0 = j["Solver"]["Fluid"]["K0"].get<double>();
            solverCompInfos.K0prime = j["Solver"]["Fluid"]["K0prime"].get<double>();
            solverCompInfos.pInfty = j["Solver"]["Fluid"]["pInfty"].get<double>();
            solverCompInfos.solverInfos = std::move(solverInfos);

            SolverCompressible solver(std::move(solverCompInfos));
            addExctractors(solver, j);
            j.clear();
            solver.solveProblem(verboseOutput);
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
