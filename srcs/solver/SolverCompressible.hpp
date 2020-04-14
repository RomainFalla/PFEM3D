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

        /**
         * \brief Apply boundary conditions to the matrix Mrho and vector Frho
         * \param invMrho A reference to the inverse mass matrix.
         * \param Frho A reference to the vector Frho.
         */
        void applyBoundaryConditionsCont(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invMrho, Eigen::VectorXd& Frho);

        /**
         * \brief Apply boundary conditions to the matrix invM and vector F
         * \param invM A reference to the inverse mass matrix.
         * \param F A reference to the vector F.
         * \param qVPrev A reference to the precedent velocity of all nodes.
         */
        void applyBoundaryConditionsMom(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F, const Eigen::VectorXd& qVPrev);

        /**
         * \brief Build the rhs of the continuity equation.
         * \param Frho A reference to the vector Fho.
         */
        void buildFrho(Eigen::VectorXd& Frho);

        /**
         * \brief Build the matrix Mrho and the vector Frho.
         * \param invMrho A reference to the inverse mass matrix.
         * \param Frho A reference to the vector Frho.
         */
        void buildMatricesCont(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invMrho, Eigen::VectorXd& Frho);

        /**
         * \brief Build the matrix M and the vector F.
         * \param invM A reference to the inverse mass matrix.
         * \param F A reference to the vector F.
         */
        void buildMatricesMom(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F);

        /**
         * \param qRho The vector of nodal density.
         * \return The vector of nodal pressure, following a Tait-Murnagham state equation.
         */
        inline Eigen::VectorXd getPFromRhoTaitMurnagham(Eigen::VectorXd qRho) const;
};

#include "SolverCompressible.inl"

#endif // SOLVERCOMPRESSIBLE_HPP_INCLUDED
