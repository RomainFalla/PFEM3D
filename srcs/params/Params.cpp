#include <iostream>
#include <fstream>

#if defined(_OPENMP)
    #include <cstdlib>
    #include <omp.h>
#endif

#include <Eigen/Core>

#include "Params.hpp"

Params::Params()
{

}

Params::~Params()
{

}

void Params::loadFromFile(std::string fileName)
{
    std::ifstream paramFile(fileName);

    if(!paramFile.is_open())
    	throw std::runtime_error("The params file cannot be read!");

    paramFile >> m_json;
    paramFile.close();

    // set the desired number of OpenMP threads
#if defined(_OPENMP)
    const char* pNumThreads = std::getenv("OMP_NUM_THREADS");

    if(pNumThreads == nullptr)
        m_threadsOMP = 1;
    else
        m_threadsOMP = std::atoi(pNumThreads);

    omp_set_num_threads(m_threadsOMP);
    Eigen::setNbThreads(m_threadsOMP);
#else
    m_threadsOMP = 1;
#endif
}
