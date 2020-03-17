#ifndef SOLVERCOMPRESSIBLE_HPP_INCLUDED
#define SOLVERCOMPRESSIBLE_HPP_INCLUDED

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../mesh/Mesh.hpp"

/**
 * \struct FluidCompressibleParams
 * \brief Incompressible fluid parameters.
 */
struct FluidCompressibleParams
{
    double rho0; /**< The fluid density (kg/m^3). */
    double mu;  /**< The fluid viscosity (Pa s). */
    double K0;
    double K0prime;
    double pInfty;
};

/**
 * \struct TimeCompressibleParams
 * \brief Time integration parameters.
 */
struct TimeCompressibleParams
{
    bool adaptDT;               /**< Should the time step be changed during the computation? */
    double currentDT;           /**< Current time step. */
    double securityCoeff;
    double maxDT;               /**< Maximu allowed time step. */
    double simuTime;            /**< Physical time which we want to simulate. */
    double simuDTToWrite;       /**< For which time interval should the program write data. */
    double nextWriteTrigger;    /**< Next physical time at which the program will write data. */
    double currentTime;
    unsigned int currentStep;
};

/**
 * \struct SolverCompressibleParams
 * \brief Compressible solver parameters.
 */
struct SolverCompressibleParams
{
    double gravity;         /**< Acceleration of the gravity (g > 0). */
    bool strongContinuity;
    FluidCompressibleParams fluid;      /**< Fluid parameters. */
    TimeCompressibleParams time;        /**< Time integration parameters. */
    std::array<bool, 6> whatToWrite {false}; /**< Which data will be written (u, v, p, ke, velocity, rho). */
    std::string writeAs;
    std::string resultsName;    /**< File name in which the results will be written. */
    std::vector<double> initialCondition;    /**< Initial condition on u, v and p (mainly to have inlet) **/
#if defined(_OPENMP)
    unsigned int nOMPThreads;
#endif
};

/**
 * \class SolverCompressible
 * \brief Represents a solver for an compressible Newtonian fluid.
 */
class SolverCompressible
{
    public:
        SolverCompressible(const nlohmann::json& j, std::string mshName, std::string resultsName);
        ~SolverCompressible();

        /**
         * \brief Display the parameters in SolverCompressibleParams structure.
         */
        void displaySolverParams() const;

        /**
         * \return A Eigen vector containing the states of the node.
         * \param beginState The first state which will be contained in q.
         * \param endState The last state which will be contained in q.
         */
        inline Eigen::VectorXd getQFromNodesStates(unsigned short beginState, unsigned short endState) const;


        /**
         * \brief Set the initial condition on u, v, p for the initial cloud of nodes.
         */
        void setInitialCondition();

        /**
         * \brief Set the states of the nodes from a Eigen vector.
         * \param q a Eigen vector containing the value of the new states. Its size is (endState-beginState+1).
         * \param beginState The first state contained in q.
         * \param endState The last state contained in q.
         */
        void setNodesStatesfromQ(const Eigen::VectorXd& q, unsigned short beginState, unsigned short endState);

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
        SolverCompressibleParams m_p;     /**< Solver parameters. */
        bool m_verboseOutput;               /**< Should the output be verbose? */

        Mesh m_mesh;                       /**< The mesh the solver is using. */

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
         * \param elementIndex The index of the element in the element list.
         * \return The gradient shape function matrix for the element in the format:
         *         [dN1dx dN2dx dN3dx 0 0 0; 0 0 0 dN1dy dN2dy dN3dy; dN1dx dN2dx dN3dx dN1dy dN2dy dN3dy]
         */
        inline Eigen::MatrixXd getB(std::size_t elementIndex) const;

        /**
         * \return The gradient shape function matrix for each gauss point in the format:
         *         [N1 N2 N3 0 0 0; 0 0 0 N1 N2 N3]
         */
        inline std::vector<Eigen::MatrixXd> getN() const;

        /**
         * \param qRho The vector of nodal density.
         * \return The vector of nodal pressure, following a Tait-Murnagham state equation.
         */
        inline Eigen::VectorXd getPFromRhoTaitMurnagham(Eigen::VectorXd qRho) const;
};

#include "SolverCompressible.inl"

#endif // SOLVERCOMPRESSIBLE_HPP_INCLUDED
