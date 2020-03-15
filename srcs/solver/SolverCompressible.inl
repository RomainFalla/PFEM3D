#include "SolverCompressible.hpp"

inline Eigen::MatrixXd SolverCompressible::getB(std::size_t elementIndex) const
{
    assert(elementIndex < m_mesh.getElementsNumber() && "elementIndex should be between 0 and size - 1 !");

    Eigen::MatrixXd B(3, 6); B.setZero();

    B(0,0) = B(2,0) = - m_mesh.getElementInvJ(elementIndex, 0, 0) - m_mesh.getElementInvJ(elementIndex, 1, 0);
    B(0,1) = B(2,1) = m_mesh.getElementInvJ(elementIndex, 0, 0);
    B(0,2) = B(2,2) = m_mesh.getElementInvJ(elementIndex, 1, 0);

    B(1,3) = B(2,3) = - m_mesh.getElementInvJ(elementIndex, 0, 1) - m_mesh.getElementInvJ(elementIndex, 1, 1);
    B(1,4) = B(2,4) = m_mesh.getElementInvJ(elementIndex, 0, 1);
    B(1,5) = B(2,5) = m_mesh.getElementInvJ(elementIndex, 1, 1);

    return B;
}

inline std::vector<Eigen::MatrixXd> SolverCompressible::getN() const
{
    std::vector<Eigen::MatrixXd> Ns;

    for(auto point: GP2Dpoints<double>)
    {
        Eigen::MatrixXd N(2, 6); N.setZero();

        N(0,0) = N(1,3) =  1 - point[0] - point[1];
        N(0,1) = N(1,4) =  point[0];
        N(0,2) = N(1,5) =  point[1];

        Ns.push_back(N);
    }

    return Ns;
}

inline Eigen::VectorXd SolverCompressible::getPFromRhoTaitMurnagham(Eigen::VectorXd qRho) const
{
    Eigen::VectorXd qP(m_mesh.getNodesNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        double rho = qRho[n];
        qP[n] = m_p.fluid.pInfty + (m_p.fluid.K0/m_p.fluid.K0prime)*(std::pow(rho/m_p.fluid.rho0, m_p.fluid.K0prime)-1);
    }

    return qP;
}

inline Eigen::VectorXd SolverCompressible::getQFromNodesStates(unsigned short beginState, unsigned short endState) const
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
