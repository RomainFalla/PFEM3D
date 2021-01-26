#include "MinMaxExtractor.hpp"

#include "../Problem.hpp"


MinMaxExtractor::MinMaxExtractor(Problem* pProblem, const std::string& outFileName, double timeBetweenWriting,
                                 unsigned short coordinate, const std::string& minMax) :
Extractor(pProblem, outFileName, timeBetweenWriting),
m_coordinate(coordinate),
m_minMax(minMax)
{
    if(m_coordinate > m_pProblem->getMesh().getDim() - 1)
        throw std::runtime_error("bad coordinate to write!");

    if(m_minMax != "min" && m_minMax != "max")
        throw std::runtime_error("unknown data to compute for MinMaxExtractor.");

    m_outFile.open(m_outFileName);
    if(!m_outFile.is_open())
        throw std::runtime_error("cannot open file to write point extractor: " + m_outFileName);
}

MinMaxExtractor::~MinMaxExtractor()
{
    m_outFile.close();
}

void MinMaxExtractor::update()
{
    if(m_pProblem->getCurrentSimTime() < m_nextWriteTrigger)
        return;

    const Mesh& mesh = m_pProblem->getMesh();

    double valueToWrite;

    if(m_minMax == "min")
    {
        valueToWrite = std::numeric_limits<double>::max();
        for(std::size_t n = 0 ; n < mesh.getNodesCount() ; ++n)
        {
            const Node& node = mesh.getNode(n);
            if(!node.isBound() && !node.isFree())
                valueToWrite = std::min(valueToWrite, node.getCoordinate(m_coordinate));
        }
    }
    else
    {
        valueToWrite = 0;
        for(std::size_t n = 0 ; n < mesh.getNodesCount() ; ++n)
        {
            const Node& node = mesh.getNode(n);
            if(!node.isBound() && !node.isFree())
                valueToWrite = std::max(valueToWrite, node.getCoordinate(m_coordinate));
        }
    }

    m_outFile << std::to_string(m_pProblem->getCurrentSimTime()) << "," << std::to_string(valueToWrite) << std::endl;

    m_nextWriteTrigger += m_timeBetweenWriting;
}
