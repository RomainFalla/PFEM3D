#pragma once
#ifndef SOLVERCOMPRESSIBLE_HPP_INCLUDED
#define SOLVERCOMPRESSIBLE_HPP_INCLUDED

#include "Solver.hpp"

/**
 * \class SolverCompressible
 * \brief Represents a solver for an compressible Newtonian fluid.
 */
class SOLVER_API SolverCompressible : public Solver
{
    public:
        SolverCompressible(const nlohmann::json& j, const std::string& mshName);
        ~SolverCompressible();

        /**
         * \brief Display the parameters in SolverCompressibleParams structure.
         */
        void displaySolverParams() const;

        /**
         * \brief Solve the Picard algorithm for one time step.
         * \return true if the algorithm converged, false otherwise.
         */
        bool solveCurrentTimeStep();

        /**
         * \brief Solve the problem for a certain set of parameters.
         */
        void solveProblem();

    private:
        double m_rho0; /**< The fluid density (kg/m^3). */
        double m_mu;  /**< The fluid viscosity (Pa s). */
        double m_K0;
        double m_K0prime;
        double m_pInfty;

        double m_securityCoeff;

        bool m_strongContinuity;

        Eigen::MatrixXd m_sumNTN;
        Eigen::MatrixXd m_ddev;

        Eigen::VectorXd m_qVPrev;           /**< The precedent speed. */
        Eigen::VectorXd m_qAccPrev;         /**< The precedent acceleration. */
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> m_invM;    /**< The mass matrix for momentum equation.. */
        Eigen::VectorXd m_F;                /**< The rhs of the momentum equation. */
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> m_invMrho; /**< The mass matrix of the continuity. */
        Eigen::VectorXd m_Frho;             /**< The rhos of the continuity equation.. */

        /**
         * \brief Apply boundary conditions to the matrix M and vector F
         */
        void applyBoundaryConditionsMom();

        /**
         * \brief Apply boundary conditions to the matrix Mrho and vector Frho
         */
        void applyBoundaryConditionsCont();

        /**
         * \brief Build the rhs of the continuity equation.
         */
        void buildFrho();

        /**
         * \brief Build the matrix M and the vector F.
         */
        void buildMatricesMom();

        /**
         * \brief Build the matrix Mrho and the vector Frho.
         */
        void buildMatricesCont();

        /**
         * \param qRho The vector of nodal density.
         * \return The vector of nodal pressure, following a Tait-Murnagham state equation.
         */
        inline Eigen::VectorXd getPFromRhoTaitMurnagham(Eigen::VectorXd qRho) const;

        /**
         * \brief Set the initial condition on u, v, p for the initial cloud of nodes.
         */
        void setInitialCondition();

};

#include "SolverCompressible.inl"

#endif // SOLVERCOMPRESSIBLE_HPP_INCLUDED
