#include "MassExtractor.hpp"

#include <vector>

#include "../Problem.hpp"


MassExtractor::MassExtractor(Problem* pProblem, const std::string& outFileName, double timeBetweenWriting) :
Extractor(pProblem, outFileName, timeBetweenWriting)
{
    m_outFile.open(m_outFileName);
    if(!m_outFile.is_open())
        throw std::runtime_error("cannot open file to write point extractor: " + m_outFileName);

    std::vector<std::string> globalData = m_pProblem->getGlobalWrittableDataName();

    if(std::find(globalData.begin(), globalData.end(), std::string("mass")) == globalData.end())
        throw std::runtime_error("the problem " + m_pProblem->getID() + " has no concept of mass !");
}

MassExtractor::~MassExtractor()
{
    m_outFile.close();
}

void MassExtractor::update()
{
    if(m_pProblem->getCurrentSimTime() < m_nextWriteTrigger)
        return;

    double valueToWrite = m_pProblem->getGlobalWrittableData(std::string("mass"));

    m_outFile << std::to_string(m_pProblem->getCurrentSimTime()) << "," << std::to_string(valueToWrite) << std::endl;

    m_nextWriteTrigger += m_timeBetweenWriting;
}
