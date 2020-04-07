#include "SolverCompressible.hpp"

inline Eigen::VectorXd SolverCompressible::getPFromRhoTaitMurnagham(Eigen::VectorXd qRho) const
{
    Eigen::VectorXd qP(m_mesh.getNodesNumber());

    #pragma omp parallel for default(shared)
    for(IndexType n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        double rho = qRho[n];
        qP[n] = m_pInfty + (m_K0/m_K0prime)*(std::pow(rho/m_rho0, m_K0prime)-1);
    }

    return qP;
}
