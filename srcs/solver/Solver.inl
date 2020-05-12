#include "Solver.hpp"

inline Eigen::MatrixXd Solver::getB(IndexType elementIndex) const noexcept
{
    assert(elementIndex < m_mesh.getElementsNumber());

    Eigen::MatrixXd B;
    if(m_mesh.getDim() == 2)
    {
        B.resize(3, 6); B.setZero();
        B(0,0) = B(2,3) = - m_mesh.getElementInvJ(elementIndex, 0, 0) - m_mesh.getElementInvJ(elementIndex, 1, 0);
        B(0,1) = B(2,4) = m_mesh.getElementInvJ(elementIndex, 0, 0);
        B(0,2) = B(2,5) = m_mesh.getElementInvJ(elementIndex, 1, 0);

        B(1,3) = B(2,0) = - m_mesh.getElementInvJ(elementIndex, 0, 1) - m_mesh.getElementInvJ(elementIndex, 1, 1);
        B(1,4) = B(2,1) = m_mesh.getElementInvJ(elementIndex, 0, 1);
        B(1,5) = B(2,2) = m_mesh.getElementInvJ(elementIndex, 1, 1);
    }
    else
    {
        B.resize(6, 12); B.setZero();
        B(0,0) = B(3,4) = B(4,8) = - m_mesh.getElementInvJ(elementIndex, 0, 0) - m_mesh.getElementInvJ(elementIndex, 1, 0) - m_mesh.getElementInvJ(elementIndex, 2, 0);
        B(0,1) = B(3,5) = B(4,9) = m_mesh.getElementInvJ(elementIndex, 0, 0);
        B(0,2) = B(3,6) = B(4,10) = m_mesh.getElementInvJ(elementIndex, 1, 0);
        B(0,3) = B(3,7) = B(4,11) = m_mesh.getElementInvJ(elementIndex, 2, 0);

        B(1,4) = B(3,0) = B(5,8) = - m_mesh.getElementInvJ(elementIndex, 0, 1) - m_mesh.getElementInvJ(elementIndex, 1, 1) - m_mesh.getElementInvJ(elementIndex, 2, 1);
        B(1,5) = B(3,1) = B(5,9) = m_mesh.getElementInvJ(elementIndex, 0, 1);
        B(1,6) = B(3,2) = B(5,10) = m_mesh.getElementInvJ(elementIndex, 1, 1);
        B(1,7) = B(3,3) = B(5,11) = m_mesh.getElementInvJ(elementIndex, 2, 1);

        B(2,8) = B(4,0) = B(5,4) = - m_mesh.getElementInvJ(elementIndex, 0, 2) - m_mesh.getElementInvJ(elementIndex, 1, 2) - m_mesh.getElementInvJ(elementIndex, 2, 2);
        B(2,9) = B(4,1) = B(5,5) = m_mesh.getElementInvJ(elementIndex, 0, 2);
        B(2,10) = B(4,2) = B(5,6) = m_mesh.getElementInvJ(elementIndex, 1, 2);
        B(2,11) = B(4,3) = B(5,7) = m_mesh.getElementInvJ(elementIndex, 2, 2);
    }

    return B;
}

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

inline Eigen::VectorXd Solver::getElementState(IndexType elm, unsigned short state) const noexcept
{
    assert(elm < m_mesh.getElementsNumber());
    assert(state < m_statesNumber);

    Eigen::VectorXd stateVec(m_mesh.getElement(elm).size());

    if(m_mesh.getDim() == 2)
    {
        stateVec << m_mesh.getNodeState(m_mesh.getElement(elm)[0], state),
                    m_mesh.getNodeState(m_mesh.getElement(elm)[1], state),
                    m_mesh.getNodeState(m_mesh.getElement(elm)[2], state);
    }
    else
    {
        stateVec << m_mesh.getNodeState(m_mesh.getElement(elm)[0], state),
                    m_mesh.getNodeState(m_mesh.getElement(elm)[1], state),
                    m_mesh.getNodeState(m_mesh.getElement(elm)[2], state),
                    m_mesh.getNodeState(m_mesh.getElement(elm)[3], state);
    }

    return stateVec;
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

inline std::vector<Eigen::MatrixXd> Solver::getN() const noexcept
{
    std::vector<Eigen::MatrixXd> Ns;

    for(unsigned short p = 0 ; p < m_mesh.getGaussPointsNumber() ; ++p)
    {
        if(m_mesh.getDim() == 2)
        {
            Eigen::MatrixXd N(2, 6); N.setZero();

            N(0,0) = N(1,3) =  1 - m_mesh.getGaussPoints(p, 0) - m_mesh.getGaussPoints(p, 1);
            N(0,1) = N(1,4) =  m_mesh.getGaussPoints(p, 0);
            N(0,2) = N(1,5) =  m_mesh.getGaussPoints(p, 1);

            Ns.push_back(std::move(N));
        }
        else
        {
            Eigen::MatrixXd N(3, 12); N.setZero();

            N(0,0) = N(1,4) = N(2,8)  =  1 - m_mesh.getGaussPoints(p, 0) - m_mesh.getGaussPoints(p, 1) - m_mesh.getGaussPoints(p, 2);
            N(0,1) = N(1,5) = N(2,9)  =  m_mesh.getGaussPoints(p, 0);
            N(0,2) = N(1,6) = N(2,10) =  m_mesh.getGaussPoints(p, 1);
            N(0,3) = N(1,7) = N(2,11) =  m_mesh.getGaussPoints(p, 2);

            Ns.push_back(std::move(N));
        }
    }

    return Ns;
}

inline Eigen::VectorXd Solver::getQFromNodesStates(unsigned short beginState, unsigned short endState) const noexcept
{
    assert(beginState < m_statesNumber && endState < m_statesNumber);

    Eigen::VectorXd q((endState - beginState + 1)*m_mesh.getNodesNumber());

    #pragma omp parallel for default(shared)
    for(IndexType n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        for (unsigned short s = beginState ; s <= endState ; ++s)
            q((s - beginState)*m_mesh.getNodesNumber() + n) = m_mesh.getNodeState(n, s);
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
    assert(q.rows() == (endState - beginState + 1)*m_mesh.getNodesNumber());

    #pragma omp parallel for default(shared)
    for(IndexType n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        for (unsigned short s = beginState ; s <= endState ; ++s)
            m_mesh.setNodeState(n, s, q((s - beginState)*m_mesh.getNodesNumber() + n));
    }
}

template <typename T, typename... Args>
void Solver::addExtractor(Args... args)
{
    static unsigned short GMSHExtractorCount = 0;
    if constexpr (std::is_same<T, PointExtractor>::value)
    {
        static_assert(sizeof...(Args) == 4, "PointExtractor needs two arguments!");
        m_pExtractor.push_back(std::make_unique<PointExtractor>(*this,
                                                                args...));
    }
    else if constexpr (std::is_same<T, MinMaxExtractor>::value)
    {
        static_assert(sizeof...(Args) == 4, "MinMaxExtractor needs two arguments ");
        m_pExtractor.push_back(std::make_unique<MinMaxExtractor>(*this,
                                                                 args...));
    }
    else if constexpr (std::is_same<T, MassExtractor>::value)
    {
        static_assert(sizeof...(Args) == 2, "MassExtractor needs 0 arguments!");
        m_pExtractor.push_back(std::make_unique<MassExtractor>(*this,
                                                               args...));
    }
    else if constexpr (std::is_same<T, GMSHExtractor>::value)
    {
        static_assert(sizeof...(Args) == 4, "GMSHExtractor needs 4 arguments!");
        if(GMSHExtractorCount == 0)
        {
            m_pExtractor.push_back(std::make_unique<GMSHExtractor>(*this,
                                                                   args...));

            GMSHExtractorCount ++;
        }
        else
            throw std::runtime_error("cannot add more than one GMSHExtractor!");
    }
    else
        static_assert(true, "Unknown extractor type ");
}
