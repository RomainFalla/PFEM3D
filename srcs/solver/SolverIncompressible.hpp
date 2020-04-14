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

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solverLU; /**< Eigen SparseLU solver. */

        /**
         * \brief Apply boundary conditions to the matrix A.
         * \param A the global sparse matrix
         * \param b the global rhs
         * \param qprev a vector containing the solution at the previous time step
         */
        void applyBoundaryConditions(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::VectorXd& qprev);

        /**
         * \brief Build the matrix A and the vector b of the Picard Algorithm.
         * \param A the global sparse matrix
         * \param b the global rhs
         * \param M the mass matrix
         * \param K the viscosity matrix
         * \param D the pressure gradient matrix
         * \param F the body force vector
         * \param tauPSPG a vector containing the value of tauPSPG for each element
         * \param qprev a vector containing the solution at the previous time step
         */
        void buildPicardSystem(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b,
                               Eigen::SparseMatrix<double>& M, Eigen::SparseMatrix<double>& K,
                               Eigen::SparseMatrix<double>& D, Eigen::VectorXd& F, const std::vector<double>& tauPSPG,
                               const Eigen::VectorXd& qprev);

        /**
         * \brief Build the coefficient tauPSPG for each element.
         * \param tauPSPG a vector which will contain the value of tauPSPG for each element
         */
        void computeTauPSPG(std::vector<double>& tauPSPG);
};

#endif // SOLVERINCOMPRESSIBLE_HPP_INCLUDED
