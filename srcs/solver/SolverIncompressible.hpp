#pragma once
#ifndef SOLVERINCOMPRESSIBLE_HPP_INCLUDED
#define SOLVERINCOMPRESSIBLE_HPP_INCLUDED

#include "Solver.hpp"

struct SolverIncompCreateInfo
{
    double rho = 0;
    double mu = 0;
    double picardRelTol = 0;
    unsigned int picardMaxIter = 0;
    double coeffDTincrease = 1;
    double coeffDTdecrease = 1;
    SolverCreateInfo solverInfos = {};
};

/**
 * \class SolverIncompressible
 * \brief Represents a solver for an incompressible Newtonian fluid.
 */
class SOLVER_API SolverIncompressible : public Solver
{
    public:
        SolverIncompressible()                                                            = delete;
        SolverIncompressible(const SolverIncompCreateInfo& solverIncompInfos);
        SolverIncompressible(const SolverIncompressible& solverIncompressible)            = delete;
        SolverIncompressible& operator=(const SolverIncompressible& solverIncompressible) = delete;
        SolverIncompressible(SolverIncompressible&& solverIncompressible)                 = delete;
        SolverIncompressible& operator=(SolverIncompressible&& solverIncompressible)      = delete;

        /**
         * \brief Display the parameters in SolverIncompressibleParams structure.
         */
        void displaySolverParams() const noexcept;

        /**
         * \brief Solve the Picard algorithm for one time step.
         * \return true if the algorithm converged, false otherwise.
         */
        bool solveCurrentTimeStep(bool verboseOutput);

        /**
         * \brief Solve the problem for a certain set of parameters.
         */
        void solveProblem(bool verboseOutput);

    private:
        double m_rho; /**< The fluid density (kg/m^3). */
        double m_mu;  /**< The fluid viscosity (Pa s). */
        double m_picardRelTol;                  /**< Relative tolerance of the algorithm. */
        unsigned int m_picardMaxIter;           /**< Maximum number of iterations for the algorithm. */
        unsigned int m_picardCurrentNumIter;    /**< Current number of Picard algorithm iterations. */
        double m_coeffDTincrease;
        double m_coeffDTdecrease;

        Eigen::MatrixXd m_MPrev;
        Eigen::VectorXd m_FPrev;
        Eigen::MatrixXd m_ddev;

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solverLU; /**< Eigen SparseLU solver. */

        /**
         * \brief Apply boundary conditions to the vector b.
         * \param b the global rhs
         * \param qprev a vector containing the solution at the previous time step
         */
        void applyBoundaryConditions(Eigen::VectorXd& b, const Eigen::VectorXd& qprev);

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
        void buildPicardSystem(Eigen::SparseMatrix<double>& A,
                               Eigen::VectorXd& b,
                               Eigen::SparseMatrix<double>& M,
                               Eigen::SparseMatrix<double>& K,
                               Eigen::SparseMatrix<double>& D,
                               Eigen::SparseMatrix<double>& C,
                               Eigen::VectorXd& F,
                               Eigen::VectorXd& H,
                               const std::vector<double>& tauPSPG, const Eigen::VectorXd& qPrev,
                               bool verboseOutput);

        /**
         * \brief Build the coefficient tauPSPG for each element.
         * \param tauPSPG a vector which will contain the value of tauPSPG for each element
         */
        void computeTauPSPG(std::vector<double>& tauPSPG);
};

#endif // SOLVERINCOMPRESSIBLE_HPP_INCLUDED
