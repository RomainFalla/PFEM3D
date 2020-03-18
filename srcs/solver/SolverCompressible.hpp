#ifndef SOLVERCOMPRESSIBLE_HPP_INCLUDED
#define SOLVERCOMPRESSIBLE_HPP_INCLUDED

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../mesh/Mesh.hpp"
#include "Solver.hpp"

/**
 * \class SolverCompressible
 * \brief Represents a solver for an compressible Newtonian fluid.
 */
class SolverCompressible : public Solver
{
    public:
        SolverCompressible(const nlohmann::json& j, const std::string& mshName, const std::string& resultsName);
        ~SolverCompressible();

        /**
         * \brief Display the parameters in SolverCompressibleParams structure.
         */
        void displaySolverParams() const;

        /**
         * \brief Set the initial condition on u, v, p for the initial cloud of nodes.
         */
        void setInitialCondition();

        /**
         * \brief Solve the Picard algorithm for one time step.
         * \return true if the algorithm converged, false otherwise.
         */
        bool solveCurrentTimeStep();

        /**
         * \brief Solve the problem for a certain set of parameters.
         */
        void solveProblem();

        /**
         * \brief Write solutions to a file.
         */
        void writeData() const;

    private:
        bool m_verboseOutput;               /**< Should the output be verbose? */

        double m_rho0; /**< The fluid density (kg/m^3). */
        double m_mu;  /**< The fluid viscosity (Pa s). */
        double m_K0;
        double m_K0prime;
        double m_pInfty;

        double m_securityCoeff;

        bool m_strongContinuity;
        std::array<bool, 6> m_whatToWrite {false}; /**< Which data will be written (u, v, p, ke, velocity, rho). */


        std::vector<Eigen::MatrixXd> m_N;     /**< The shape functions matrices for eah gauss point (wihout *detJ). */
        Eigen::MatrixXd m_sumNTN;
        Eigen::VectorXd m_m;
        Eigen::MatrixXd m_ddev;

        Eigen::VectorXd m_qVPrev;           /**< The precedent speed. */
        Eigen::VectorXd m_qAccPrev;         /**< The precedent acceleration. */
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> m_invM;    /**< The mass matrix for momentum equation.. */
        Eigen::VectorXd m_F;                /**< The rhs of the momentum equation. */
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> m_invMrho; /**< The mass matrix of the continuity. */
        Eigen::VectorXd m_Frho;             /**< The rhos of the continuity equation.. */

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solverLU; /**< Eigen SparseLU solver. */

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
};

#include "SolverCompressible.inl"

#endif // SOLVERCOMPRESSIBLE_HPP_INCLUDED
