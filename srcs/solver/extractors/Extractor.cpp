#include "Extractor.hpp"

Extractor::Extractor(const Solver& solver, const std::string& outFileName, double timeBetweenWriting) :
m_solver(solver),
m_outFileName(outFileName),
m_timeBetweenWriting(timeBetweenWriting),
m_nextWriteTrigger(0)
{
    if(m_timeBetweenWriting <= 0)
        throw std::runtime_error("the interval between write should be strictly greater than 0!");
}

Extractor::~Extractor()
{
}

void Extractor::update()
{
    throw std::runtime_error("Extractor base class is useless ^^!");
}
