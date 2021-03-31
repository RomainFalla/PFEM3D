#include "PicardAlgo.hpp"


PicardAlgo::PicardAlgo(std::function<void(Eigen::SparseMatrix<double>& /** A **/,
                                          Eigen::VectorXd& /** b **/,
                                          const Eigen::VectorXd& /** qPrev **/)> buildAb,
                       std::function<void(Eigen::VectorXd& /** b **/,
                                          const Eigen::VectorXd& /** qPrev **/)> applyBC,
                       std::function<void(const Eigen::VectorXd& /** q **/)> exetuteTask,
                       std::function<double(const Eigen::VectorXd& /** qIter **/,
                                            const Eigen::VectorXd& /** qIterPrev **/)> computeRes,
                       unsigned int maxIter, double minRes):
NonLinearAlgo(buildAb, applyBC, exetuteTask, computeRes),
m_maxIter(maxIter),
m_minRes(minRes)
{

}

PicardAlgo::~PicardAlgo()
{

}

void PicardAlgo::displayParams()
{
    std::cout << " * Maximum residual: " << m_minRes << "\n"
              << " * Maximum iteration count: " << m_maxIter << std::endl;
}

bool PicardAlgo::solve(Mesh* pMesh, const Eigen::VectorXd& qPrev, bool verboseOutput)
{
    Eigen::SparseMatrix<double> A(qPrev.rows(), qPrev.rows());
    Eigen::VectorXd b(qPrev.rows());

    pMesh->saveNodesList();

    unsigned int iterCount = 0;
    Eigen::VectorXd qIter(qPrev.rows()); qIter.setZero();
    Eigen::VectorXd qIterPrev(qPrev.rows());
    double res = std::numeric_limits<double>::max();

    while(res > m_minRes)
    {
        if(verboseOutput)
        {
            std::cout << " - Picard algorithm (mesh position) - iteration ("
                      << iterCount << ")" << std::endl;
        }

        qIterPrev = qIter;

        m_buildAb(A, b, qPrev);
        A.makeCompressed();
        m_applyBC(b, qPrev);
        m_solver.compute(A);

        if(m_solver.info() == Eigen::Success)
        {
            qIter = m_solver.solve(b);
        }
        else
        {
            if(verboseOutput)
                std::cout << "\t * The Eigen::SparseLU solver failed to factorize the A matrix!" << std::endl;
            pMesh->restoreNodesList();
            return false;
        }

        m_executeTask(qIter);

        iterCount++;

        if(iterCount == 1)
        {
            res = std::numeric_limits<double>::max();
        }
        else
        {
            res = m_computeRes(qIter, qIterPrev);
        }

        if(verboseOutput)
            std::cout << "\t * Relative 2-norm of q: " << res << " vs "
                      << m_minRes << std::endl;


        if(iterCount > m_maxIter)
        {
            if(verboseOutput)
            {
                std::cout << "\t * Iteration count " << iterCount
                      << " greater than maximum: " << m_maxIter << std::endl;
            }

            pMesh->restoreNodesList();
            return false;
        }

        if(std::isnan(res))
        {
            if(verboseOutput)
                std::cout << "\t * NaN residual" << std::endl;

            pMesh->restoreNodesList();
            return false;
        }
    }

    return true;
}
