#include "NonLinearAlgo.hpp"

#include <stdexcept>

NonLinearAlgo::NonLinearAlgo(std::function<void(Eigen::SparseMatrix<double>& /** A **/,
                                                Eigen::VectorXd& /** b **/,
                                                const Eigen::VectorXd& /** qPrev **/)> buildAb,
                             std::function<void(Eigen::VectorXd& /** b **/,
                                                const Eigen::VectorXd& /** qPrev **/)> applyBC,
                             std::function<void(const Eigen::VectorXd& /** q **/)> executeTask,
                             std::function<double(const Eigen::VectorXd& /** qIter **/,
                                                  const Eigen::VectorXd& /** qIterPrev **/)> computeRes):
m_buildAb(buildAb),
m_applyBC(applyBC),
m_executeTask(executeTask),
m_computeRes(computeRes)
{

}

NonLinearAlgo::~NonLinearAlgo()
{

}

void NonLinearAlgo::displayParams()
{
    throw std::runtime_error("Unimplemented function by the child class -> NonLinearAlgo::displayParams()");
}

bool NonLinearAlgo::solve(Mesh* pMesh, const Eigen::VectorXd& qPrev, bool verboseOutput)
{
    (void)pMesh;
    (void)qPrev;
    (void)verboseOutput;
    throw std::runtime_error("Unimplemented function by the child class -> NonLinearAlgo::solve(pMesh, qPrev)");
}

