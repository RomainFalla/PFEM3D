#include "SolverBoussinesqWC.hpp"

inline Eigen::VectorXd SolverBoussinesqWC::getPFromRhoTaitMurnagham(Eigen::VectorXd qRho) const
{
    Eigen::VectorXd qP(m_mesh.getNodesCount());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesCount() ; ++n)
    {
        double rho = qRho[n];
        qP[n] = m_pInfty + (m_K0/m_K0prime)*(std::pow(rho/m_rho0, m_K0prime)-1);
    }

    return qP;
}
