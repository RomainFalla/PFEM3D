#include "Solver.hpp"

inline Eigen::VectorXd Solver::getQFromNodesStates(unsigned short beginState, unsigned short endState) const
{
    Eigen::VectorXd q((endState - beginState + 1)*m_mesh.getNodesNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        for (unsigned short s = beginState ; s <= endState ; ++s)
            q((s - beginState)*m_mesh.getNodesNumber() + n) = m_mesh.getNodeState(n, s);
    }

    return q;
}
