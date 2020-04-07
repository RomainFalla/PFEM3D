#pragma once
#ifndef SOLVERINCOMPRESSIBLE_HPP_INCLUDED
#define SOLVERINCOMPRESSIBLE_HPP_INCLUDED

#include "Solver.hpp"

/**
 * \class SolverIncompressible
 * \brief Represents a solver for an incompressible Newtonian fluid.
 */
class SOLVER_API SolverIncompressible : public Solver
{
    public:
        SolverIncompressible(const nlohmann::json& j, const std::string& mshName);
        ~SolverIncompressible();

        /**
         * \brief Display the parameters in SolverIncompressibleParams structure.
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
        double m_rho; /**< The fluid density (kg/m^3). */
        double m_mu;  /**< The fluid viscosity (Pa s). */
        double m_picardRelTol;                  /**< Relative tolerance of the algorithm. */
        unsigned int m_picardMaxIter;           /**< Maximum number of iterations for the algorithm. */
        unsigned int m_picardCurrentNumIter;    /**< Current number of Picard algorithm iterations. */
        double m_coeffDTincrease;
        double m_coeffDTdecrease;

        Eigen::MatrixXd m_sumNTN;
        Eigen::MatrixXd m_ddev;

        Eigen::VectorXd m_qprev;            /**< The precedent solution */
        Eigen::SparseMatrix<double> m_A;    /**< The matrix A representing the problem: [M/dt+K -D^T; C/dt-D L]. */
        Eigen::VectorXd m_b;                /**< The vector b representing the problem: [M/dt*qprev + F; H]. */
        Eigen::SparseMatrix<double> m_M;    /**< The mass matrix. */
        Eigen::SparseMatrix<double> m_K;    /**< The viscosity matrix. */
        Eigen::SparseMatrix<double> m_D;    /**< The pressure matrix. */
        Eigen::VectorXd m_F;                /**< The volume force vector. */
        std::vector<double> m_tauPSPG;      /**< tau_PSPG parameters for each element. */

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solverLU; /**< Eigen SparseLU solver. */

        /**
         * \brief Apply boundary conditions to the matrix A
         */
        void applyBoundaryConditions();

        /**
         * \brief Build the matrix A and the vector b of the Picard Algorithm.
         */
        void buildPicardSystem();

        /**
         * \brief Build the matrix A and the vector b of the Picard Algorithm.
         */
        void computeTauPSPG();

        /**
         * \brief Set the initial condition on u, v, p for the initial cloud of nodes.
         */
        void setInitialCondition();
};

#endif // SOLVERINCOMPRESSIBLE_HPP_INCLUDED
