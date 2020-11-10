#include "Solver.hpp"

inline double Solver::getCurrentDT() const noexcept
{
    return m_currentDT;
}

inline unsigned int Solver::getCurrentStep() const noexcept
{
    return m_currentStep;
}

inline double Solver::getCurrentTime() const noexcept
{
    return m_currentTime;
}

inline double Solver::getEndTime() const noexcept
{
    return m_endTime;
}

inline double Solver::getMaxDT() const noexcept
{
    return m_maxDT;
}

inline const Mesh& Solver::getMesh() const noexcept
{
    return m_mesh;
}

inline Eigen::VectorXd Solver::getQFromNodesStates(unsigned short beginState, unsigned short endState) const noexcept
{
    assert(beginState < m_statesNumber && endState < m_statesNumber);

    std::size_t nodeCount = m_mesh.getNodesCount();

    Eigen::VectorXd q((endState - beginState + 1)*nodeCount);

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < nodeCount ; ++n)
    {
        const Node& node = m_mesh.getNode(n);

        for (unsigned short s = beginState ; s <= endState ; ++s)
            q((s - beginState)*nodeCount + n) = node.getState(s);
    }

    return q;
}

inline SOLVER_TYPE Solver::getSolverType() const noexcept
{
    return m_solverType;
}

inline unsigned short Solver::getStatesNumber() const noexcept
{
    return m_statesNumber;
}

inline bool Solver::isDTAdaptable() const noexcept
{
    return m_adaptDT;
}

inline void Solver::setCurrentDT(double dt)
{
    if(!m_adaptDT)
        throw std::runtime_error("Cannot change dt of a solver which does not authorize it");

    if(dt <=0)
        throw std::runtime_error("The new time interval should be strictly positive!");

    m_currentDT = dt;
}

inline void Solver::setNodesStatesfromQ(const Eigen::VectorXd& q, unsigned short beginState, unsigned short endState) noexcept
{
    assert(beginState < m_statesNumber && endState < m_statesNumber);
    assert(static_cast<std::size_t>(q.rows()) == (endState - beginState + 1)*m_mesh.getNodesCount());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesCount() ; ++n)
    {
        for (unsigned short s = beginState ; s <= endState ; ++s)
            m_mesh.setNodeState(n, s, q((s - beginState)*m_mesh.getNodesCount() + n));
    }
}
