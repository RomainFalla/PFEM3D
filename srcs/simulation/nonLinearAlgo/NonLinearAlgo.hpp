#pragma once
#ifndef NONLINEARALGO_HPP_INCLUDED
#define NONLINEARALGO_HPP_INCLUDED

#include <functional>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "../../mesh/Mesh.hpp"

/**
 * \class NonLinearAlgo
 * \brief Represents a non-linear algorithm to solve a Ax = b system of equations.
 */
class NonLinearAlgo
{
    public:
        NonLinearAlgo(std::function<void(Eigen::SparseMatrix<double>& /** A **/,
                                         Eigen::VectorXd& /** b **/,
                                         const Eigen::VectorXd& /** qPrev **/)> buildAb,
                      std::function<void(Eigen::VectorXd& /** b **/,
                                         const Eigen::VectorXd& /** qPrev **/)> applyBC,
                      std::function<void(const Eigen::VectorXd& /** qIter **/)> executeTask,
                      std::function<double(const Eigen::VectorXd& /** qIter **/,
                                           const Eigen::VectorXd& /** qIterPrev **/)> computeRes);

        virtual ~NonLinearAlgo();

        virtual void displayParams();
        virtual bool solve(Mesh* pMesh, const Eigen::VectorXd& qPrev, bool verboseOutput);

    protected:
        std::function<void(Eigen::SparseMatrix<double>&,
                           Eigen::VectorXd&,
                           const Eigen::VectorXd&)> m_buildAb;

        std::function<void(Eigen::VectorXd&,
                           const Eigen::VectorXd&)> m_applyBC;

        std::function<void(const Eigen::VectorXd&)> m_executeTask;

        std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&)> m_computeRes;
};

#endif // NONLINEARALGO_HPP_INCLUDED
