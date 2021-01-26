#pragma once
#ifndef PICARDALGO_HPP_INCLUDED
#define PICARDALGO_HPP_INCLUDED

#include <functional>
#include <Eigen/SparseLU>

#include "NonLinearAlgo.hpp"

/**
 * \class PicardAlgo
 * \brief Represents a Picard non-linear algorithm to solve a Ax = b system of equations.
 */
class PicardAlgo : public NonLinearAlgo
{
    public:
        PicardAlgo(std::function<void(Eigen::SparseMatrix<double>& /** A **/,
                                      Eigen::VectorXd& /** b **/,
                                      const Eigen::VectorXd& /** qPrev **/)> buildAb,
                   std::function<void(Eigen::VectorXd& /** b **/,
                                      const Eigen::VectorXd& /** qPrev **/)> applyBC,
                   std::function<void(const Eigen::VectorXd& /** qIter**/)> executeTask,
                   std::function<double(const Eigen::VectorXd& /** qIter **/,
                                        const Eigen::VectorXd& /** qIterPrev **/)> computeRes,
                   unsigned int maxIter, double minRes);

        ~PicardAlgo() override;

        void displayParams() override;
        bool solve(Mesh* pMesh, const Eigen::VectorXd& qPrev, bool verboseOutput) override;
    private:
        unsigned int m_maxIter;
        double m_minRes;

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solver;
};

#endif // PICARDALGO_HPP_INCLUDED
