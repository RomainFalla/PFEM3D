#include "Problem.hpp"

inline const Mesh& Problem::getMesh() const noexcept
{
    return *m_pMesh;
}

inline const Solver& Problem::getSolver() const noexcept
{
    return *m_pSolver;
}

inline double Problem::getCurrentSimTime() const noexcept
{
    return m_time;
}

inline double Problem::getCurrentSimStep() const noexcept
{
    return m_step;
}

inline unsigned int Problem::getThreadCount() const noexcept
{
    return m_nThreads;
}

inline bool Problem::isOutputVerbose() const noexcept
{
    return m_verboseOutput;
}


