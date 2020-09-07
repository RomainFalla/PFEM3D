#include "GMSHExtractor.hpp"

#include <utility>
#include <gmsh.h>

#include "../Solver.hpp"


GMSHExtractor::GMSHExtractor(const Solver& solver, const std::string& outFileName, double timeBetweenWriting,
                             const std::vector<std::string>& whatToWrite, std::string writeAs) :
Extractor(solver, outFileName, timeBetweenWriting),
m_writeAs(std::move(writeAs))
{
    if(!(m_writeAs == "Nodes" || m_writeAs == "Elements" || m_writeAs == "NodesElements"))
        throw std::runtime_error("unknown data type for results writing " + m_writeAs);

    if(solver.getSolverType() == SOLVER_TYPE::Incompressible_PSPG)
    {
        if(solver.getMesh().getDim() == 2)
            m_whatCanBeWritten = {"u", "v", "p", "ke", "velocity"};
        else
            m_whatCanBeWritten = {"u", "v", "w", "p", "ke", "velocity"};
    }
    else if(solver.getSolverType() == SOLVER_TYPE::WeaklyCompressible)
    {
        if(solver.getMesh().getDim() == 2)
            m_whatCanBeWritten = {"u", "v", "p", "rho", "ax", "ay", "ke", "velocity"};
        else
            m_whatCanBeWritten = {"u", "v", "w", "p", "rho", "ax", "ay", "az", "ke", "velocity"};
    }
    else throw std::runtime_error("unknown solver type!");

    m_whatToWrite.resize(m_whatCanBeWritten.size());

    for(unsigned short i = 0 ; i < whatToWrite.size() ; ++i)
    {
        bool found = false;
        for(unsigned short j = 0 ; j < m_whatCanBeWritten.size() ; ++j)
        {
            if(whatToWrite[i] == m_whatCanBeWritten[j])
            {
                m_whatToWrite[j] = true;
                found = true;
                break;
            }
        }
        if(!found)
            throw std::runtime_error("unknown quantity to write " + whatToWrite[i]);
    }

    gmsh::initialize();
#ifndef NDEBUG
    gmsh::option::setNumber("General.Terminal", 1);
#else
    gmsh::option::setNumber("General.Terminal", 0);
#endif // DEBUG
}

GMSHExtractor::~GMSHExtractor()
{
    gmsh::finalize();
}

void GMSHExtractor::update()
{
    if(m_solver.getCurrentTime() < m_nextWriteTrigger)
        return;

    const Mesh& mesh = m_solver.getMesh();

    gmsh::model::add("theModel");
    gmsh::model::setCurrent("theModel");
    if(m_writeAs == "Nodes" || m_writeAs == "NodesElements")
        gmsh::model::addDiscreteEntity(0, 1);
    if(m_writeAs == "Elements" || m_writeAs == "NodesElements")
        gmsh::model::addDiscreteEntity(mesh.getDim(), 2);

    for(unsigned short i = 0 ; i < m_whatToWrite.size() ; ++i)
    {
        if(m_whatToWrite[i])
        gmsh::view::add(m_whatCanBeWritten[i], i + 1);
    }

    std::vector<std::size_t> nodesTags(mesh.getNodesNumber());
    std::vector<double> nodesCoord(3*mesh.getNodesNumber());
    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < mesh.getNodesNumber() ; ++n)
    {
        nodesTags[n] = n + 1;
        nodesCoord[3*n] = mesh.getNodePosition(n, 0);
        nodesCoord[3*n + 1] = mesh.getNodePosition(n, 1);
        if(mesh.getDim() == 2)
            nodesCoord[3*n + 2] = 0;
        else
            nodesCoord[3*n + 2] = mesh.getNodePosition(n, 2);
    }

    gmsh::model::mesh::addNodes(0, 1, nodesTags, nodesCoord);

    if(m_writeAs == "Elements" || m_writeAs == "NodesElements")
    {
        std::vector<std::size_t> elementTags(mesh.getElementsNumber());
        std::vector<std::size_t> nodesTagsPerElement((mesh.getDim()+1)*mesh.getElementsNumber());
        #pragma omp parallel for default(shared)
        for(std::size_t elm = 0 ; elm < mesh.getElementsNumber() ; ++elm)
        {
            elementTags[elm] = elm + 1;
            for(std::size_t n = 0 ; n < mesh.getElement(elm).size() ; ++n)
                nodesTagsPerElement[(mesh.getDim() + 1)*elm + n] = mesh.getElement(elm)[n] + 1;
        }

        if(mesh.getDim() == 2)
            gmsh::model::mesh::addElementsByType(2, 2, elementTags, nodesTagsPerElement); //Triangle 2
        else
            gmsh::model::mesh::addElementsByType(2, 4, elementTags, nodesTagsPerElement); //Tetrahedron 4
    }

    if(m_writeAs == "Nodes" || m_writeAs == "NodesElements")
        gmsh::model::mesh::addElementsByType(1, 15, nodesTags, nodesTags);

    std::vector<std::vector<std::vector<double>>> dataNodesStates;
    std::vector<std::vector<double>> dataKe;
    std::vector<std::vector<double>> dataVelocity;

    dataNodesStates.resize(m_solver.getStatesNumber());

    for(unsigned short i = 0 ; i < m_solver.getStatesNumber() ; ++i)
    {
        if(m_whatToWrite[i])
            dataNodesStates[i].resize(mesh.getNodesNumber());
    }

    if(m_whatToWrite[m_solver.getStatesNumber()])
        dataKe.resize(mesh.getNodesNumber());
    if(m_whatToWrite[m_solver.getStatesNumber() + 1])
        dataVelocity.resize(mesh.getNodesNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < mesh.getNodesNumber() ; ++n)
    {
        for(unsigned short i = 0 ; i < m_solver.getStatesNumber() ; ++i)
        {
            if(m_whatToWrite[i])
            {
                const std::vector<double> state{mesh.getNodeState(n, i)};
                dataNodesStates[i][n] = state;
            }
        }
        if(m_whatToWrite[m_solver.getStatesNumber()])
        {
            double ke = 0;
            for(unsigned short d = 0 ; d < mesh.getDim() ; ++d)
                ke += mesh.getNodeState(n, d)*mesh.getNodeState(n, d);

            ke = 0.5*std::sqrt(ke);

            const std::vector<double> keVector{ke};
            dataKe[n] = keVector;
        }
        if(m_whatToWrite[m_solver.getStatesNumber() + 1])
        {
            if(mesh.getDim() == 2)
            {
                const std::vector<double> velocity{mesh.getNodeState(n, 0), mesh.getNodeState(n, 1), 0};
                dataVelocity[n] = velocity;
            }
            else
            {
                const std::vector<double> velocity{mesh.getNodeState(n, 0), mesh.getNodeState(n, 1), mesh.getNodeState(n, 2)};
                dataVelocity[n] = velocity;
            }
        }
    }

    for(unsigned short i = 0 ; i < m_solver.getStatesNumber() ; ++i)
    {
        if(m_whatToWrite[i])
            gmsh::view::addModelData(i + 1, m_solver.getCurrentStep(), "theModel", "NodeData", nodesTags, dataNodesStates[i], m_solver.getCurrentTime(), 1);
    }

    if(m_whatToWrite[m_solver.getStatesNumber()])
        gmsh::view::addModelData(m_solver.getStatesNumber() + 1, m_solver.getCurrentStep(), "theModel", "NodeData", nodesTags, dataKe, m_solver.getCurrentTime(), 1);
    if(m_whatToWrite[m_solver.getStatesNumber() + 1])
        gmsh::view::addModelData(m_solver.getStatesNumber() + 2, m_solver.getCurrentStep(), "theModel", "NodeData", nodesTags, dataVelocity, m_solver.getCurrentTime(), 3);

    const std::string baseName = m_outFileName.substr(0, m_outFileName.find(".msh"));

    for(unsigned short i = 0; i < m_whatToWrite.size(); ++i)
    {
        if(m_whatToWrite[i] == true)
            gmsh::view::write(i + 1, baseName + "_" + std::to_string(m_solver.getCurrentTime()) + ".msh" , true);
    }

    gmsh::model::remove();

    m_nextWriteTrigger += m_timeBetweenWriting;
}
