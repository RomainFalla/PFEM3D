#include "Solver.hpp"

#include <iostream>

#include <gmsh.h>


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

    m_N = getN();

    if(m_mesh.getMeshDim() == 2)
    {
        m_m.resize(3);
        m_m << 1, 1, 0;

        m_bodyForces.resize(2);
        m_bodyForces << 0 , -m_gravity;
    }
    else if(m_mesh.getMeshDim() == 3)
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

void Solver::writeData() const
{
    assert(m_whatCanBeWriten.size() == m_statesNumber + 2);
    assert(m_whatCanBeWriten.size() == m_whatToWrite.size());

    gmsh::model::add("theModel");
    gmsh::model::setCurrent("theModel");
    if(m_writeAs == "Nodes" || m_writeAs == "NodesElements")
        gmsh::model::addDiscreteEntity(0, 1);
    if(m_writeAs == "Elements" || m_writeAs == "NodesElements")
        gmsh::model::addDiscreteEntity(m_mesh.getMeshDim(), 2);

    for(unsigned short i = 0 ; i < m_whatToWrite.size() ; ++i)
    {
        if(m_whatToWrite[i])
        gmsh::view::add(m_whatCanBeWriten[i], i + 1);
    }

    std::vector<std::size_t> nodesTags(m_mesh.getNodesNumber());
    std::vector<double> nodesCoord(3*m_mesh.getNodesNumber());
    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        nodesTags[n] = n + 1;
        nodesCoord[3*n] = m_mesh.getNodePosition(n, 0);
        nodesCoord[3*n + 1] = m_mesh.getNodePosition(n, 1);
        if(m_mesh.getMeshDim() == 2)
            nodesCoord[3*n + 2] = 0;
        else if(m_mesh.getMeshDim() == 3)
            nodesCoord[3*n + 2] = m_mesh.getNodePosition(n, 2);
    }

    gmsh::model::mesh::addNodes(0, 1, nodesTags, nodesCoord);

    if(m_writeAs == "Elements" || m_writeAs == "NodesElements")
    {
        std::vector<std::size_t> elementTags(m_mesh.getElementsNumber());
        std::vector<std::size_t> nodesTagsPerElement((m_mesh.getMeshDim()+1)*m_mesh.getElementsNumber());
        #pragma omp parallel for default(shared)
        for(std::size_t elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
        {
            elementTags[elm] = elm + 1;
            for(std::size_t n = 0 ; n < m_mesh.getElement(elm).size() ; ++n)
                nodesTagsPerElement[(m_mesh.getMeshDim() + 1)*elm + n] = m_mesh.getElement(elm)[n] + 1;
        }

        if(m_mesh.getMeshDim() == 2)
            gmsh::model::mesh::addElementsByType(2, 2, elementTags, nodesTagsPerElement); //Triangle 2
        else if(m_mesh.getMeshDim() == 3)
            gmsh::model::mesh::addElementsByType(2, 4, elementTags, nodesTagsPerElement); //Tetrahedron 4
    }

    if(m_writeAs == "Nodes" || m_writeAs == "NodesElements")
        gmsh::model::mesh::addElementsByType(1, 15, nodesTags, nodesTags);

    std::vector<std::vector<std::vector<double>>> dataNodesStates;
    std::vector<std::vector<double>> dataKe;
    std::vector<std::vector<double>> dataVelocity;

    dataNodesStates.resize(m_statesNumber);

    for(unsigned short i = 0 ; i < m_statesNumber ; ++i)
    {
        if(m_whatToWrite[i])
            dataNodesStates[i].resize(m_mesh.getNodesNumber());
    }

    if(m_whatToWrite[m_statesNumber])
        dataKe.resize(m_mesh.getNodesNumber());
    if(m_whatToWrite[m_statesNumber + 1])
        dataVelocity.resize(m_mesh.getNodesNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        for(unsigned short i = 0 ; i < m_statesNumber ; ++i)
        {
            if(m_whatToWrite[i])
            {
                const std::vector<double> state{m_mesh.getNodeState(n, i)};
                dataNodesStates[i][n] = state;
            }
        }
        if(m_whatToWrite[m_statesNumber])
        {
            double ke = 0;
            for(unsigned short d = 0 ; d < m_mesh.getMeshDim() ; ++d)
                ke += m_mesh.getNodeState(n, d)*m_mesh.getNodeState(n, d);
            ke /= static_cast<double>(m_mesh.getMeshDim());
            const std::vector<double> keVector{ke};
            dataKe[n] = keVector;
        }
        if(m_whatToWrite[m_statesNumber + 1])
        {
            if(m_mesh.getMeshDim() == 2)
            {
                const std::vector<double> velocity{m_mesh.getNodeState(n, 0), m_mesh.getNodeState(n, 1), 0};
                dataVelocity[n] = velocity;
            }
            else if(m_mesh.getMeshDim() == 3)
            {
                const std::vector<double> velocity{m_mesh.getNodeState(n, 0), m_mesh.getNodeState(n, 1), m_mesh.getNodeState(n, 2)};
                dataVelocity[n] = velocity;
            }
        }
    }

    for(unsigned short i = 0 ; i < m_statesNumber ; ++i)
    {
        if(m_whatToWrite[i])
            gmsh::view::addModelData(i + 1, m_currentStep, "theModel", "NodeData", nodesTags, dataNodesStates[i], m_currentTime, 1);
    }

    if(m_whatToWrite[m_statesNumber])
        gmsh::view::addModelData(m_statesNumber + 1, m_currentStep, "theModel", "NodeData", nodesTags, dataKe, m_currentTime, 1);
    if(m_whatToWrite[m_statesNumber + 1])
        gmsh::view::addModelData(m_statesNumber + 2, m_currentStep, "theModel", "NodeData", nodesTags, dataVelocity, m_currentTime, 3);

    const std::string baseName = m_resultsName.substr(0, m_resultsName.find(".msh"));

    for(unsigned short i = 0; i < m_whatToWrite.size(); ++i)
    {
        if(m_whatToWrite[i] == true)
            gmsh::view::write(i + 1, baseName + "_" + std::to_string(m_currentTime) + ".msh" , true);
    }

    gmsh::model::remove();
}
