#include "Solver.hpp"

#include <iostream>
#include <type_traits>


Solver::Solver(const SolverCreateInfo& solverInfos) :
m_gravity(solverInfos.gravity),
m_strongPAtFS(solverInfos.strongPAtFS),
m_adaptDT(solverInfos.adaptDT),
m_currentDT(solverInfos.initialDT),
m_currentStep(0), m_currentTime(0),
m_endTime(solverInfos.endTime),
m_maxDT(solverInfos.maxDT),
m_mesh(solverInfos.meshInfos),
m_hasGMSHExtractor(false)
{
    m_solverType = SOLVER_TYPE::Undefined;

    m_lua.open_libraries(sol::lib::base, sol::lib::math);

    // set the desired number of OpenMP threads
#if defined(_OPENMP)
    const char* pNumThreads = std::getenv("OMP_NUM_THREADS");

    if(pNumThreads == nullptr)
        m_numOMPThreads = 1;
    else
        m_numOMPThreads = std::atoi(pNumThreads);

    omp_set_num_threads(m_numOMPThreads);
    Eigen::setNbThreads(m_numOMPThreads);
#else
     m_numOMPThreads = 1;
#endif

    if(m_currentDT <= 0)
        throw std::runtime_error("the initial time step should be strictly greater than 0!");

    if(m_endTime <= 0)
        throw std::runtime_error("the total time to simulate should be strictly greater than 0!");

    if(m_maxDT <= 0)
        throw std::runtime_error("the maximum time interval should be strictly greater than 0!");

    if(solverInfos.IBCfile.empty())
        throw std::runtime_error("no IBC lua file provided!");

    m_lua.script_file(solverInfos.IBCfile);
    m_lua["g"] = m_gravity;

    m_N = getN();

    if(m_mesh.getDim() == 2)
    {
        m_m.resize(3);
        m_m << 1, 1, 0;

        m_bodyForces.resize(2);
        m_bodyForces << 0 , -m_gravity;
    }
    else
    {
        m_m.resize(6);
        m_m << 1, 1, 1, 0, 0, 0;

        m_bodyForces.resize(3);
        m_bodyForces << 0 , 0, -m_gravity;
    }
}

void Solver::setInitialCondition()
{
    assert(m_mesh.getNodesNumber() != 0);

    //With lua, not thread safe !
    for(IndexType n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        std::pair<std::vector<double>, bool> result;
        result = m_lua[std::string("init") + m_mesh.getNodeType(n)](m_mesh.getNodePosition(n)).get<std::pair<std::vector<double>, bool>>();

        if(result.first.size() != m_statesNumber)
            throw std::runtime_error("Your initial condition does not set the right number of state!");

        for(unsigned short i = 0 ; i < m_statesNumber ; ++i)
        {
            m_mesh.setNodeState(n, i, result.first[i]);
            m_mesh.setNodeIsDirichlet(n, result.second);
        }

    }
}

void Solver::addPointExtractor(const std::string& outFileName, double timeBetweenWriting,
                               unsigned short stateToWrite, const std::vector<std::vector<double>>& points)
{
    m_pExtractor.push_back(std::make_unique<PointExtractor>(*this,
                                                            outFileName,
                                                            timeBetweenWriting,
                                                            stateToWrite,
                                                            points));
}

void Solver::addMinMaxExtractor(const std::string& outFileName, double timeBetweenWriting,
                                unsigned short coordinate, const std::string& minMax)
{
    m_pExtractor.push_back(std::make_unique<MinMaxExtractor>(*this,
                                                             outFileName,
                                                             timeBetweenWriting,
                                                             coordinate,
                                                             minMax));
}

void Solver::addMassExtractor(const std::string& outFileName, double timeBetweenWriting)
{
    m_pExtractor.push_back(std::make_unique<MassExtractor>(*this,
                                                           outFileName,
                                                           timeBetweenWriting));
}

void Solver::addGMSHExtractor(const std::string& outFileName, double timeBetweenWriting,
                              const std::vector<std::string>& whatToWrite, std::string writeAs)
{
    if(m_hasGMSHExtractor)
    {
        std::cerr << "There is already a GMSH extractor set" << std::endl;
        return;
    }

    m_pExtractor.push_back(std::make_unique<GMSHExtractor>(*this,
                                                           outFileName,
                                                           timeBetweenWriting,
                                                           whatToWrite,
                                                           writeAs));

    m_hasGMSHExtractor = true;
}
