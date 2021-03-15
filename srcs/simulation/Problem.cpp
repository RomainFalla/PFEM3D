#include "Problem.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <omp.h>
#include <Eigen/Core>

#include "Solver.hpp"
#include "extractors/Extractors.hpp"
#include "../mesh/Mesh.hpp"

Problem::Problem(const std::string& luaFilePath):
m_time(0),
m_step(0)
{
    std::cout   << "================================================================"
                << "\n"
                << "                        LOADING PROBLEM                         "
                << "\n"
                << "================================================================"
                << std::endl;

    const char* pNumThreads = std::getenv("OMP_NUM_THREADS");

    if(pNumThreads == nullptr)
        m_nThreads = 1;
    else
        m_nThreads = std::atoi(pNumThreads);

    omp_set_num_threads(m_nThreads);
    Eigen::setNbThreads(m_nThreads);

    m_solState.resize(m_nThreads);
    m_problemParams.resize(m_nThreads);

    for(std::size_t i = 0 ; i < m_nThreads ; ++i)
    {
        m_solState[i].open_libraries(sol::lib::base, sol::lib::math, sol::lib::table);
        m_solState[i].script_file(luaFilePath);

        m_problemParams[i] = SolTable("Problem", m_solState[i]);
    }

    m_id = m_problemParams[0].checkAndGet<std::string>("id");

    SolTable mesh = SolTable("Mesh", m_problemParams[0]);

    MeshCreateInfo createInfo = {};

    createInfo.hchar = mesh.checkAndGet<double>("hchar");
    createInfo.alpha = mesh.checkAndGet<double>("alpha");
    createInfo.gamma = mesh.checkAndGet<double>("gamma");
    createInfo.omega = mesh.checkAndGet<double>("omega");
    createInfo.useMeshRefinement = mesh.checkAndGet<bool>("useMeshRefinement");
    createInfo.boundingBox = mesh.checkAndGet<std::vector<double>>("boundingBox");
    createInfo.mshFile = mesh.checkAndGet<std::string>("mshFile");

    if (createInfo.useMeshRefinement) 
    {
        SolTable refinement = SolTable("Refinement", mesh);

        createInfo.refInfo = {};

        createInfo.refInfo.alphaRatio = refinement.checkAndGet<double>("alphaRatio");
        createInfo.refInfo.gamma2 = refinement.checkAndGet<double>("gamma2");
        createInfo.refInfo.minTargetMeshSize = refinement.checkAndGet<double>("minTargetMeshSize");
        createInfo.refInfo.maxTargetMeshSize = refinement.checkAndGet<double>("maxTargetMeshSize");
        createInfo.refInfo.maxProgressionFactor = refinement.checkAndGet<double>("maxProgressionFactor");

    }


    std::cout << "Loading the mesh" << std::flush;
    m_pMesh = std::make_unique<Mesh>(createInfo);
    std::cout << "\rLoading the mesh\t\tok" << std::endl;

    m_maxTime = m_problemParams[0].checkAndGet<double>("simulationTime");
    m_verboseOutput = m_problemParams[0].checkAndGet<bool>("verboseOutput");
}

Problem::~Problem()
{

}

void Problem::displayParams() const
{
    throw std::runtime_error("Unimplemented function by the child class -> Problem::displayParams()");
}

std::string Problem::getID() const noexcept
{
    return m_id;
}

std::vector<std::string> Problem::getWrittableDataName() const
{
    throw std::runtime_error("Unimplemented function by the child class -> Problem::getWrittableDataName()");
}

std::vector<double> Problem::getWrittableData(const std::string& /** name **/, std::size_t /** nodeIndex **/) const
{
    throw std::runtime_error("Unimplemented function by the child class -> Problem::getWrittableData(name, nodeIndex)");
}

std::vector<std::string> Problem::getGlobalWrittableDataName() const
{
    throw std::runtime_error("Unimplemented function by the child class -> Problem::getGlobalWrittableDataName()");
}

double Problem::getGlobalWrittableData(const std::string& /** name **/) const
{
    throw std::runtime_error("Unimplemented function by the child class -> Problem::getGlobalWrittableData(name)");
}

std::vector<std::string> Problem::getMeshWrittableDataName() const
{
    return {"normals", "curvatures", "debug"};
}

std::vector<double> Problem::getMeshWrittableData(const std::string& name, std::size_t nodeIndex) const
{
    if(name == "normals")
    {
        if(!m_pMesh->isNormalCurvComputed())
            return {0, 0, 0};

        if(m_pMesh->getNode(nodeIndex).isOnFreeSurface() || m_pMesh->getNode(nodeIndex).isBound())
        {
            std::array<double, 3> normal = m_pMesh->getBoundFSNormal(nodeIndex);
            return {normal[0], normal[1], normal[2]};
        }
        else
            return {0, 0, 0};
    }
    else if(name == "curvatures")
    {
        if(!m_pMesh->isNormalCurvComputed())
            return {0};

        if(m_pMesh->getNode(nodeIndex).isOnFreeSurface())
        {
            double curvature = m_pMesh->getFreeSurfaceCurvature(nodeIndex);
            return {curvature};
        }
        else
            return {0};
    }
    else if(name == "debug")
    {
        const Node& node = m_pMesh->getNode(nodeIndex);

        if(node.isOnFreeSurface())
            return {-2};
        else if(node.isFree())
            return {-1};
        else if(node.isBound() || node.isContact())
            return {1};
        else if(node.isFixed())
            return {2};
        else
            return {0};
    }
    else
        throw std::runtime_error("The mesh data " + name + " cannot be extracted from the mesh");
}

void Problem::addExtractors()
{
    std::cout << "Loading extractors" << std::flush;

    SolTable extractors = SolTable("Extractors", m_problemParams[0]);

    extractors.for_each([this](sol::object /*key*/, sol::object value){
        SolTable extractor = SolTable(value);
        std::string kind = extractor.checkAndGet<std::string>("kind");

        std::string outFileName = extractor.checkAndGet<std::string>("outputFile");
        double timeBetweenWriting = extractor.checkAndGet<double>("timeBetweenWriting");

        if(kind == "GMSH")
        {
            std::vector<std::string> whatToWrite = extractor.checkAndGet<std::vector<std::string>>("whatToWrite");
            std::string writeAs = extractor.checkAndGet<std::string>("writeAs");

            m_pExtractors.push_back(std::make_unique<GMSHExtractor>(this,
                                                                    outFileName,
                                                                    timeBetweenWriting,
                                                                    whatToWrite,
                                                                    writeAs));
        }
        else if(kind == "Point")
        {
            std::string whatToWrite = extractor.checkAndGet<std::string>("whatToWrite");
            std::vector<std::vector<double>> points = extractor.checkAndGet<std::vector<std::vector<double>>>("points");

            m_pExtractors.push_back(std::make_unique<PointExtractor>(this,
                                                                     outFileName,
                                                                     timeBetweenWriting,
                                                                     whatToWrite,
                                                                     points));
        }
        else if(kind == "Mass")
        {
            m_pExtractors.push_back(std::make_unique<MassExtractor>(this,
                                                                    outFileName,
                                                                    timeBetweenWriting));
        }
        else if(kind == "MinMax")
        {
            unsigned short coordinate = extractor.checkAndGet<unsigned short>("coordinate");
            std::string minMax = extractor.checkAndGet<std::string>("minMax");

            m_pExtractors.push_back(std::make_unique<MinMaxExtractor>(this,
                                                                      outFileName,
                                                                      timeBetweenWriting,
                                                                      coordinate,
                                                                      minMax));
        }
        else
            throw std::runtime_error("Unknown extractor kind " + kind + "!");
    });

    std::cout << "\rLoading extractors\t\tok" << std::endl;
}

void Problem::setInitialCondition()
{
    std::cout << "Setting initial conditions" << std::flush;

    SolTable initialCond = SolTable("IC", m_problemParams[0]);
    initialCond.checkCall("initStates", m_pMesh->getNode(0).getPosition());

    assert(m_pMesh->getNodesCount() != 0);

    //With lua, not thread safe !
    for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);

        std::vector<double> result;
        result = initialCond.call<std::vector<double>>("initStates", node.getPosition());

        if(result.size() != m_statesNumber)
            throw std::runtime_error("Your initial condition does not set the right number of state: " +
                                     std::to_string(result.size()) + " vs " + std::to_string(m_statesNumber) + "!");


        if(node.isBound())
        {
            bool isFixed = initialCond.checkAndGet<bool>(m_pMesh->getNodeType(n) + std::string("Fixed"));
            m_pMesh->setNodeIsFixed(n, isFixed);

            bool ok = initialCond.checkCallNoThrow("init" + m_pMesh->getNodeType(n) + "States", node.getPosition());

            if(ok)
                result = initialCond.call<std::vector<double>>("init" + m_pMesh->getNodeType(n) + "States", node.getPosition());
        }

        for(unsigned short i = 0 ; i < m_statesNumber ; ++i)
            m_pMesh->setNodeState(n, i, result[i]);
    }

    std::cout << "\rSetting initial conditions\tok" << std::endl;
}

void Problem::updateTime(double timeStep)
{
    m_time += timeStep;
    m_step++;
}

void Problem::simulate()
{
    std::cout   << std::string(20 + m_id.size(), '=')
                << std::endl
                << " EXECUTING " << m_id << " problem"
                << std::endl
                << std::string(20 + m_id.size(), '=')
                << std::endl;

    displayParams();

    std::cout << std::string(40, '-') << std::endl;

    for(auto& pExtractor : m_pExtractors)
    {
        pExtractor->update();
    }

    if (m_pMesh->useMeshRefinement())
    {
        m_pMesh->updateBoundingElementSize();
        m_pMesh->updateLocalElementSize();
    }

    while(m_time < m_maxTime)
    {
        if(m_verboseOutput)
        {
            std::cout << "----------------------------------------------------------------" << std::endl;
            std::cout << "Solving time step: " << m_time + m_pSolver->getTimeStep()
                      << "/" << m_maxTime << " s, dt = " << m_pSolver->getTimeStep() << std::endl;
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
        else
        {
            if(m_step%100 == 0 || m_time + m_pSolver->getTimeStep() >= m_maxTime)
            {
                std::cout << std::fixed << std::setprecision(3);
                std::cout << "\r" << "Solving time step: " << m_time + m_pSolver->getTimeStep()
                          << "/" << m_maxTime << " s, dt = ";
                std::cout << std::scientific;
                std::cout << m_pSolver->getTimeStep() << " s" << std::flush;
            }
        }

        bool ok = m_pSolver->solveOneTimeStep();

        if(ok)
        {
            for(auto& pExtractor : m_pExtractors)
            {
                pExtractor->update();
            }
            if (m_pMesh->useMeshRefinement())
                m_pMesh->updateLocalElementSize(); // The local element size is updated at each iteration
        }

        m_pSolver->computeNextDT();
    }

    std::cout << std::endl;
}
