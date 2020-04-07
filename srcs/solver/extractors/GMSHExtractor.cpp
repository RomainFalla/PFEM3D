#include "GMSHExtractor.hpp"

#include <utility>
#include <gmsh.h>


GMSHExtractor::GMSHExtractor(const std::string& outFileName, double timeBetweenWriting,
                             const std::vector<std::string>& whatToWrite,
                             std::vector<std::string> whatCanBeWritten, std::string writeAs,
                             unsigned short statesNumber) :
Extractor(outFileName, timeBetweenWriting), m_whatCanBeWritten(std::move(whatCanBeWritten)), m_writeAs(std::move(writeAs)),
m_statesNumber(statesNumber)
{
    if(!(m_writeAs == "Nodes" || m_writeAs == "Elements" || m_writeAs == "NodesElements"))
    {
        std::string errotText = std::string("unknown data type for results writing ") + m_writeAs;
        throw std::runtime_error(errotText);
    }

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
        {
            std::string errorText = std::string("unknown quantity to write ") + whatToWrite[i];
            throw std::runtime_error(errorText);
        }
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

void GMSHExtractor::update(const Mesh& mesh, double currentTime, unsigned int currentStep)
{
    if(currentTime < m_nextWriteTrigger)
        return;

    gmsh::model::add("theModel");
    gmsh::model::setCurrent("theModel");
    if(m_writeAs == "Nodes" || m_writeAs == "NodesElements")
        gmsh::model::addDiscreteEntity(0, 1);
    if(m_writeAs == "Elements" || m_writeAs == "NodesElements")
        gmsh::model::addDiscreteEntity(mesh.getMeshDim(), 2);

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
        if(mesh.getMeshDim() == 2)
            nodesCoord[3*n + 2] = 0;
        else if(mesh.getMeshDim() == 3)
            nodesCoord[3*n + 2] = mesh.getNodePosition(n, 2);
    }

    gmsh::model::mesh::addNodes(0, 1, nodesTags, nodesCoord);

    if(m_writeAs == "Elements" || m_writeAs == "NodesElements")
    {
        std::vector<std::size_t> elementTags(mesh.getElementsNumber());
        std::vector<std::size_t> nodesTagsPerElement((mesh.getMeshDim()+1)*mesh.getElementsNumber());
        #pragma omp parallel for default(shared)
        for(std::size_t elm = 0 ; elm < mesh.getElementsNumber() ; ++elm)
        {
            elementTags[elm] = elm + 1;
            for(std::size_t n = 0 ; n < mesh.getElement(elm).size() ; ++n)
                nodesTagsPerElement[(mesh.getMeshDim() + 1)*elm + n] = mesh.getElement(elm)[n] + 1;
        }

        if(mesh.getMeshDim() == 2)
            gmsh::model::mesh::addElementsByType(2, 2, elementTags, nodesTagsPerElement); //Triangle 2
        else if(mesh.getMeshDim() == 3)
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
            dataNodesStates[i].resize(mesh.getNodesNumber());
    }

    if(m_whatToWrite[m_statesNumber])
        dataKe.resize(mesh.getNodesNumber());
    if(m_whatToWrite[m_statesNumber + 1])
        dataVelocity.resize(mesh.getNodesNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < mesh.getNodesNumber() ; ++n)
    {
        for(unsigned short i = 0 ; i < m_statesNumber ; ++i)
        {
            if(m_whatToWrite[i])
            {
                const std::vector<double> state{mesh.getNodeState(n, i)};
                dataNodesStates[i][n] = state;
            }
        }
        if(m_whatToWrite[m_statesNumber])
        {
            double ke = 0;
            for(unsigned short d = 0 ; d < mesh.getMeshDim() ; ++d)
                ke += mesh.getNodeState(n, d)*mesh.getNodeState(n, d);

            ke = 0.5*std::sqrt(ke);

            const std::vector<double> keVector{ke};
            dataKe[n] = keVector;
        }
        if(m_whatToWrite[m_statesNumber + 1])
        {
            if(mesh.getMeshDim() == 2)
            {
                const std::vector<double> velocity{mesh.getNodeState(n, 0), mesh.getNodeState(n, 1), 0};
                dataVelocity[n] = velocity;
            }
            else if(mesh.getMeshDim() == 3)
            {
                const std::vector<double> velocity{mesh.getNodeState(n, 0), mesh.getNodeState(n, 1), mesh.getNodeState(n, 2)};
                dataVelocity[n] = velocity;
            }
        }
    }

    for(unsigned short i = 0 ; i < m_statesNumber ; ++i)
    {
        if(m_whatToWrite[i])
            gmsh::view::addModelData(i + 1, currentStep, "theModel", "NodeData", nodesTags, dataNodesStates[i], currentTime, 1);
    }

    if(m_whatToWrite[m_statesNumber])
        gmsh::view::addModelData(m_statesNumber + 1, currentStep, "theModel", "NodeData", nodesTags, dataKe, currentTime, 1);
    if(m_whatToWrite[m_statesNumber + 1])
        gmsh::view::addModelData(m_statesNumber + 2, currentStep, "theModel", "NodeData", nodesTags, dataVelocity, currentTime, 3);

    const std::string baseName = m_outFileName.substr(0, m_outFileName.find(".msh"));

    for(unsigned short i = 0; i < m_whatToWrite.size(); ++i)
    {
        if(m_whatToWrite[i] == true)
            gmsh::view::write(i + 1, baseName + "_" + std::to_string(currentTime) + ".msh" , true);
    }

    gmsh::model::remove();

    m_nextWriteTrigger += m_timeBetweenWriting;
}
