#include "MinMaxExtractor.hpp"

#include <iostream>
#include <Eigen/Dense>


MinMaxExtractor::MinMaxExtractor(const std::string& outFileName, double timeBetweenWriting,
                                 unsigned short coordinate, const std::string& minMax, unsigned short meshDim) :
Extractor(outFileName, timeBetweenWriting), m_coordinate(coordinate), m_minMax(minMax)
{
    if(m_coordinate > meshDim - 1)
        throw std::runtime_error("bad coordinate to write!");

    if(m_minMax != "min" && m_minMax != "max")
        throw std::runtime_error("unknown data to compute for MinMaxExtractor.");

    m_outFile.open(m_outFileName);
    if(!m_outFile.is_open())
    {
        std::string errorText = std::string("cannot open file to write point extractor: ") + m_outFileName;
        throw std::runtime_error(errorText);
    }

}

MinMaxExtractor::~MinMaxExtractor()
{
    m_outFile.close();
}

void MinMaxExtractor::update(const Mesh& mesh, double currentTime, unsigned int currentStep)
{
    if(currentTime >= m_nextWriteTrigger)
    {
        double valueToWrite;

        if(m_minMax == "min")
        {
            valueToWrite = std::numeric_limits<double>::max();
            for(IndexType n = 0 ; n < mesh.getNodesNumber() ; ++n)
            {
                if(!mesh.isNodeBound(n) && !mesh.isNodeBound(n))
                    valueToWrite = std::min(valueToWrite, mesh.getNodePosition(n, m_coordinate));
            }
        }
        else
        {
            valueToWrite = 0;
            for(IndexType n = 0 ; n < mesh.getNodesNumber() ; ++n)
            {
                if(!mesh.isNodeBound(n) && !mesh.isNodeBound(n))
                    valueToWrite = std::max(valueToWrite, mesh.getNodePosition(n, m_coordinate));
            }
        }

        m_outFile << std::to_string(currentTime) << "," << std::to_string(valueToWrite) << std::endl;

        m_nextWriteTrigger += m_timeBetweenWriting;
    }
}
