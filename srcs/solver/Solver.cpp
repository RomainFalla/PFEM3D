#include "Solver.hpp"

Solver::Solver(const nlohmann::json& j, const std::string& mshName, const std::string& resultsName) :
m_mesh(j)
{
    m_mesh.loadFromFile(mshName);

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

    m_resultsName = resultsName;

    m_gravity                   = j["Solver"]["gravity"].get<double>();

    m_initialCondition          = j["Solver"]["initialCondition"].get<std::vector<double>>();
    m_writeAs                   = j["Solver"]["writeAs"].get<std::string>();

    if(!(m_writeAs == "Nodes" || m_writeAs == "Elements" || m_writeAs == "NodesElements"))
        throw std::runtime_error("Unknown data type for results writing!");

    m_adaptDT                   = j["Solver"]["Time"]["adaptDT"].get<bool>();
    m_currentDT                 = j["Solver"]["Time"]["initialDT"].get<double>();

    if(m_currentDT <= 0)
        throw std::runtime_error("The initial time step should be strictly greater than 0!");

    m_currentStep               = 0;
    m_currentTime               = 0;
    m_endTime                   = j["Solver"]["Time"]["endTime"].get<double>();

    if(m_endTime <= 0)
        throw std::runtime_error("The total time to simulate should be strictly greater than 0!");

    m_maxDT                     = j["Solver"]["Time"]["maxDT"].get<double>();

    if(m_maxDT <= 0)
        throw std::runtime_error("The maximum time interval should be strictly greater than 0!");

    m_timeBetweenWriting        = j["Solver"]["Time"]["timeBetweenWriting"].get<double>();

    if(m_timeBetweenWriting <= 0)
        throw std::runtime_error("The interval between write should be strictly greater than 0!");

    m_nextWriteTrigger   = m_timeBetweenWriting;
}

Solver::~Solver()
{

}
