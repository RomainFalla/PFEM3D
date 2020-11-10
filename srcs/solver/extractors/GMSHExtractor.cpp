#include "GMSHExtractor.hpp"

#include <utility>
#include <gmsh.h>

#include "../Solver.hpp"

uint16_t GMSHExtractor::m_initialized;

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
    else if(solver.getSolverType() == SOLVER_TYPE::Conduction)
    {
            m_whatCanBeWritten = {"T"};
    }
    else if(solver.getSolverType() == SOLVER_TYPE::Boussinesq)
    {
        if(solver.getMesh().getDim() == 2)
            m_whatCanBeWritten = {"u", "v", "p", "T", "ke", "velocity"};
        else
            m_whatCanBeWritten = {"u", "v", "w", "p", "T", "ke", "velocity"};
    }
    else if(solver.getSolverType() == SOLVER_TYPE::BoussinesqWC)
    {
        if(solver.getMesh().getDim() == 2)
            m_whatCanBeWritten = {"u", "v", "p", "rho",  "ax", "ay", "T", "ke", "velocity"};
        else
            m_whatCanBeWritten = {"u", "v", "w", "p", "rho", "ax", "ay", "az", "T", "ke", "velocity"};
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

    if(m_initialized == 0)
    {
        gmsh::initialize();
#ifndef NDEBUG
        gmsh::option::setNumber("General.Terminal", 1);
#else
        gmsh::option::setNumber("General.Terminal", 0);
#endif // DEBUG
    }

    m_initialized++;
}

GMSHExtractor::~GMSHExtractor()
{
    if(m_initialized == 1)
    {
        gmsh::finalize();
    }

    m_initialized--;
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

    std::vector<std::size_t> nodesTags(mesh.getNodesCount());
    std::vector<double> nodesCoord(3*mesh.getNodesCount());
    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < mesh.getNodesCount() ; ++n)
    {
        const Node& node = mesh.getNode(n);
        nodesTags[n] = n + 1;
        nodesCoord[3*n] = node.getPosition(0);
        nodesCoord[3*n + 1] = node.getPosition(1);
        if(mesh.getDim() == 2)
            nodesCoord[3*n + 2] = 0;
        else
            nodesCoord[3*n + 2] = node.getPosition(2);
    }

    gmsh::model::mesh::addNodes(0, 1, nodesTags, nodesCoord);

    if(m_writeAs == "Elements" || m_writeAs == "NodesElements")
    {
        std::vector<std::size_t> elementTags(mesh.getElementsCount());
        std::vector<std::size_t> nodesTagsPerElement((mesh.getDim()+1)*mesh.getElementsCount());
        #pragma omp parallel for default(shared)
        for(std::size_t elm = 0 ; elm < mesh.getElementsCount() ; ++elm)
        {
            const Element& element = mesh.getElement(elm);

            elementTags[elm] = elm + 1;
            for(std::size_t n = 0 ; n < element.getNodeCount() ; ++n)
                nodesTagsPerElement[(mesh.getDim() + 1)*elm + n] = element.getNodeIndex(n) + 1;
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
            dataNodesStates[i].resize(mesh.getNodesCount());
    }

    if(m_whatToWrite[m_solver.getStatesNumber()])
        dataKe.resize(mesh.getNodesCount());
    if(m_whatToWrite[m_solver.getStatesNumber() + 1])
        dataVelocity.resize(mesh.getNodesCount());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < mesh.getNodesCount() ; ++n)
    {
        const Node& node = mesh.getNode(n);

        for(unsigned short i = 0 ; i < m_solver.getStatesNumber() ; ++i)
        {
            if(m_whatToWrite[i])
            {
                const std::vector<double> state{node.getState(i)};
                dataNodesStates[i][n] = state;
            }
        }
        if(m_whatToWrite[m_solver.getStatesNumber()])
        {
            double ke = 0;
            for(unsigned short d = 0 ; d < mesh.getDim() ; ++d)
                ke += node.getState(d)*node.getState(d);

            ke = 0.5*std::sqrt(ke);

            const std::vector<double> keVector{ke};
            dataKe[n] = keVector;
        }
        if(m_whatToWrite[m_solver.getStatesNumber() + 1])
        {
            if(mesh.getDim() == 2)
            {
                const std::vector<double> velocity{node.getState(0), node.getState(1), 0};
                dataVelocity[n] = velocity;
            }
            else
            {
                if(node.isOnFreeSurface())
                {
                    std::array<double, 3> normal = mesh.getFreeSurfaceNormal(n);
                    const std::vector<double> velocity{normal[0], normal[1], normal[2]};
                    dataVelocity[n] = velocity;
                }
                else
                {
                    const std::vector<double> velocity{0, 0, 0};
                    dataVelocity[n] = velocity;
                }
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
