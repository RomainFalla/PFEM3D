#include "Solver.hpp"

#include <iostream>


Solver::Solver(const nlohmann::json& j, const std::string& mshName) :
m_mesh(j)
{
    m_solverType = Undefined;
    m_mesh.loadFromFile(mshName);

    m_lua.open_libraries(sol::lib::base, sol::lib::math);

    m_verboseOutput             = j["verboseOutput"].get<bool>();

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

    m_gravity                   = j["Solver"]["gravity"].get<double>();

    m_adaptDT                   = j["Solver"]["Time"]["adaptDT"].get<bool>();
    m_currentDT                 = j["Solver"]["Time"]["initialDT"].get<double>();

    if(m_currentDT <= 0)
        throw std::runtime_error("the initial time step should be strictly greater than 0!");

    m_currentStep               = 0;
    m_currentTime               = 0;
    m_endTime                   = j["Solver"]["Time"]["endTime"].get<double>();

    if(m_endTime <= 0)
        throw std::runtime_error("the total time to simulate should be strictly greater than 0!");

    m_maxDT                     = j["Solver"]["Time"]["maxDT"].get<double>();

    if(m_maxDT <= 0)
        throw std::runtime_error("the maximum time interval should be strictly greater than 0!");

    std::string IBConditions = j["Solver"]["IBCs"].get<std::string>();
    m_lua.script_file(IBConditions);

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

Solver::~Solver()
{

}

void Solver::setInitialCondition()
{
    assert(m_mesh.getNodesNumber() != 0);

    //With lua, not thread safe !
    for(IndexType n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        std::pair<std::vector<double>, bool> result;
        result = m_lua[std::string("init") + m_mesh.getNodeType(n)](m_mesh.getNodePosition(n)).get<std::pair<std::vector<double>, bool>>();

        for(unsigned short i = 0 ; i < m_statesNumber ; ++i)
        {
            m_mesh.setNodeState(n, i, result.first[i]);
            m_mesh.setNodeIsDirichlet(n, result.second);
        }

    }
}
