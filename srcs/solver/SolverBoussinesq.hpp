#pragma once
#ifndef SOLVERBOUSSINESQ_HPP_INCLUDED
#define SOLVERBOUSSINESQ_HPP_INCLUDED

#include "Solver.hpp"

struct SolverBoussinesqCreateInfo
{
    double rho = 0;
    double mu = 0;
    double cv = 0;
    double k = 0;
    double alpha = 0;
    double Tr = 0;
    double gamma = 0;
    double picardRelTol = 0;
    unsigned int picardMaxIter = 0;
    double coeffDTincrease = 1;
    double coeffDTdecrease = 1;
    SolverCreateInfo solverInfos = {};
};

/**
 * \class SolverBoussinesq
 * \brief Represents a solver for an incompressible Newtonian fluid with boussinesq approximation.
 */
class SOLVER_API SolverBoussinesq : public Solver
{
    public:
        SolverBoussinesq()                                                                = delete;
        SolverBoussinesq(const SolverBoussinesqCreateInfo& solverBoussinesqIncomp);
        SolverBoussinesq(const SolverBoussinesq& solverBoussinesqIncomp)            = delete;
        SolverBoussinesq& operator=(const SolverBoussinesq& solverBoussinesqIncomp) = delete;
        SolverBoussinesq(SolverBoussinesq&& solverBoussinesqIncomp)                 = delete;
        SolverBoussinesq& operator=(SolverBoussinesq&& solverBoussinesqIncomp)      = delete;

        /**
         * \brief Display the parameters of the problem.
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
        double m_gamma; /**< The fluid surface tension (N/m). */
        double m_k; /**< The fluid thermal conductivity (W/(mK)). */
        double m_cv; /**< The fluid heat capacity. */
        double m_alpha;
        double m_Tr;
        double m_picardRelTol;                  /**< Relative tolerance of the algorithm. */
        unsigned int m_picardMaxIter;           /**< Maximum number of iterations for the algorithm. */
        unsigned int m_picardCurrentNumIter;    /**< Current number of Picard algorithm iterations. */
        double m_coeffDTincrease;
        double m_coeffDTdecrease;

        enum class BC
        {
            T = 1,
            Q = 2,
            V = 3,
            TandV = 4,
            QandV = 5,
            None = -1
        };

        Eigen::MatrixXd m_MPrev;
        Eigen::MatrixXd m_MPrevT;
        Eigen::MatrixXd m_ddev;

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solverLU; /**< Eigen SparseLU solver. */

        /**
         * \brief Apply boundary conditions to the vector b (heat equation).
         * \param b the global rhs
         * \param qprev a vector containing the solution at the previous time step
         */
        void applyBoundaryConditionsMom(Eigen::VectorXd& b, const Eigen::VectorXd& qprev);

        /**
         * \brief Build the matrix A and the vector b of the Picard Algorithm for the momentum.
         * \param A the global sparse matrix
         * \param b the global rhs
         * \param qprev a vector containing the solution at the previous time step
         */
        void buildPicardSystemMom(Eigen::SparseMatrix<double>& A,
                               Eigen::VectorXd& b,
                               const Eigen::VectorXd& qPrev);

        /**
         * \brief Apply boundary conditions to the vector b.
         * \param b the global rhs
         * \param qprev a vector containing the solution at the previous time step
         */
        void applyBoundaryConditionsHeat(Eigen::VectorXd& b, const Eigen::VectorXd& qprev);

        /**
         * \brief Build the matrix A and the vector b for the heat equation.
         * \param A the global sparse matrix
         * \param b the global rhs
         * \param qprev a vector containing the solution at the previous time step
         */
        void buildAHeat(Eigen::SparseMatrix<double>& A,
                        Eigen::VectorXd& b,
                        const Eigen::VectorXd& qPrev);

        /**
         * \brief Build the coefficient tauPSPG for each element.
         * \param tauPSPG a vector which will contain the value of tauPSPG for each element
         */
        double computeTauPSPG(const Element& element);
};

#endif // SOLVERBOUSSINESQ_HPP_INCLUDED
