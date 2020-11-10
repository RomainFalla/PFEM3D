#pragma once
#ifndef SolverConduction_HPP_INCLUDED
#define SolverConduction_HPP_INCLUDED

#include "Solver.hpp"

struct SolverConductionCreateInfo
{
    double rho = 0;
    double cv = 0;
    double k = 0;
    SolverCreateInfo solverInfos = {};
};

/**
 * \class SolverConduction
 * \brief Represents a solver for an incompressible Newtonian fluid.
 */
class SOLVER_API SolverConduction : public Solver
{
    public:
        SolverConduction()                                                    = delete;
        SolverConduction(const SolverConductionCreateInfo& solverCondInfos);
        SolverConduction(const SolverConduction& solverConduction)            = delete;
        SolverConduction& operator=(const SolverConduction& solverConduction) = delete;
        SolverConduction(SolverConduction&& solverConduction)                 = delete;
        SolverConduction& operator=(SolverConduction&& solverConduction)      = delete;

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
        double m_k;  /**< The fluid viscosity (Pa s). */
        double m_cv;  /**< The fluid viscosity (Pa s). */

        enum class BC
        {
            T = 1,
            Q = 2,
            None = -1
        };

        Eigen::MatrixXd m_MPrev;

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
         * \param qprev a vector containing the solution at the previous time step
         */
        void buildPicardSystem(Eigen::SparseMatrix<double>& A,
                               Eigen::VectorXd& b,
                               const Eigen::VectorXd& qPrev);
};

#endif // SolverConduction_HPP_INCLUDED
