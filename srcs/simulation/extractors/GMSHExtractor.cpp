#include "GMSHExtractor.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>
#include <gmsh.h>

#include "../Problem.hpp"

unsigned int GMSHExtractor::m_initialized;

GMSHExtractor::GMSHExtractor(Problem* pProblem, const std::string& outFileName, double timeBetweenWriting,
                             const std::vector<std::string>& whatToWrite, std::string writeAs) :
Extractor(pProblem, outFileName, timeBetweenWriting),
m_writeAs(std::move(writeAs))
{
    if(!(m_writeAs == "Nodes" || m_writeAs == "Elements" || m_writeAs == "NodesElements"))
        throw std::runtime_error("unknown data type for results writing " + m_writeAs);

    std::vector<std::string> writtableData = m_pProblem->getWrittableDataName();
    std::vector<std::string> meshWrittableData = m_pProblem->getMeshWrittableDataName();

    for(std::string dataToWrite: whatToWrite)
    {
        if(std::find(writtableData.begin(), writtableData.end(), dataToWrite) != writtableData.end())
            m_whatToWrite.push_back(dataToWrite);
        else if(std::find(meshWrittableData.begin(), meshWrittableData.end(), dataToWrite) != meshWrittableData.end())
            m_meshToWrite.push_back(dataToWrite);
        else
            throw std::runtime_error("the problem " + m_pProblem->getID() + " cannot write data named " + dataToWrite + " using GMSHExtractor!");
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
    if(m_pProblem->getCurrentSimTime() < m_nextWriteTrigger)
        return;

    const Mesh& mesh = m_pProblem->getMesh();

    gmsh::model::add("theModel");
    gmsh::model::setCurrent("theModel");
    if(m_writeAs == "Nodes" || m_writeAs == "NodesElements")
        gmsh::model::addDiscreteEntity(0, 1);
    if(m_writeAs == "Elements" || m_writeAs == "NodesElements")
        gmsh::model::addDiscreteEntity(mesh.getDim(), 2);

    std::vector<std::size_t> nodesTags(mesh.getNodesCount());
    std::vector<double> nodesCoord(3*mesh.getNodesCount());
    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < mesh.getNodesCount() ; ++n)
    {
        const Node& node = mesh.getNode(n);
        nodesTags[n] = n + 1;
        nodesCoord[3*n] = node.getCoordinate(0);
        nodesCoord[3*n + 1] = node.getCoordinate(1);
        if(mesh.getDim() == 2)
            nodesCoord[3*n + 2] = 0;
        else
            nodesCoord[3*n + 2] = node.getCoordinate(2);
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
            for(std::size_t n = 0 ; n < mesh.getNodesPerElm() ; ++n)
                nodesTagsPerElement[(mesh.getDim() + 1)*elm + n] = element.getNodeIndex(n) + 1;
        }

        if(mesh.getDim() == 2)
            gmsh::model::mesh::addElementsByType(2, 2, elementTags, nodesTagsPerElement); //Triangle 2
        else
            gmsh::model::mesh::addElementsByType(2, 4, elementTags, nodesTagsPerElement); //Tetrahedron 4
    }

    if(m_writeAs == "Nodes" || m_writeAs == "NodesElements")
        gmsh::model::mesh::addElementsByType(1, 15, nodesTags, nodesTags);

    std::map<std::string, std::vector<std::vector<double>>> dataToWritePerNodePerState;

    std::size_t counter = 1;

    for(std::string toWrite : m_whatToWrite)
    {
        gmsh::view::add(toWrite, counter);

        dataToWritePerNodePerState[toWrite] = std::vector<std::vector<double>>(mesh.getNodesCount());

        counter++;
    }

    for(std::string toWrite : m_meshToWrite)
    {
        gmsh::view::add(toWrite, counter);

        dataToWritePerNodePerState[toWrite] = std::vector<std::vector<double>>(mesh.getNodesCount());

        counter++;
    }

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < mesh.getNodesCount() ; ++n)
    {
        for(std::string toWrite : m_whatToWrite)
        {
            std::vector<double> dataToWrite = m_pProblem->getWrittableData(toWrite, n);
            dataToWritePerNodePerState[toWrite][n] = dataToWrite;
        }

        for(std::string toWrite : m_meshToWrite)
        {
            std::vector<double> dataToWrite = m_pProblem->getMeshWrittableData(toWrite, n);
            dataToWritePerNodePerState[toWrite][n] = dataToWrite;
        }
    }

    const std::string baseName = m_outFileName.substr(0, m_outFileName.find(".msh"));

    counter = 1;
    for(std::string toWrite : m_whatToWrite)
    {
        gmsh::view::addModelData(counter, m_pProblem->getCurrentSimStep(), "theModel", "NodeData", nodesTags,
                                 dataToWritePerNodePerState[toWrite], m_pProblem->getCurrentSimTime(),
                                 dataToWritePerNodePerState[toWrite][0].size());

        gmsh::view::write(counter, baseName + "_" + std::to_string(m_pProblem->getCurrentSimTime()) + ".msh" , true);

        counter++;
    }

    for(std::string toWrite : m_meshToWrite)
    {
        gmsh::view::addModelData(counter, m_pProblem->getCurrentSimStep(), "theModel", "NodeData", nodesTags,
                                 dataToWritePerNodePerState[toWrite], m_pProblem->getCurrentSimTime(),
                                 dataToWritePerNodePerState[toWrite][0].size());

        gmsh::view::write(counter, baseName + "_" + std::to_string(m_pProblem->getCurrentSimTime()) + ".msh" , true);

        counter++;
    }

    gmsh::model::remove();

    m_nextWriteTrigger += m_timeBetweenWriting;
}
