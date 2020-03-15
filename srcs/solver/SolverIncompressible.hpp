#ifndef SOLVERINCOMPRESSIBLE_HPP_INCLUDED
#define SOLVERINCOMPRESSIBLE_HPP_INCLUDED

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../mesh/Mesh.hpp"

/**
 * \struct FluidIncompressibleParams
 * \brief Incompressible fluid parameters.
 */
struct FluidIncompressibleParams
{
    double rho; /**< The fluid density (kg/m^3). */
    double mu;  /**< The fluid viscosity (Pa s). */
};

/**
 * \struct PicardParams
 * \brief Picard algorithm parameters.
 */
struct PicardParams
{
    double relTol;                  /**< Relative tolerance of the algorithm. */
    unsigned int maxIter;           /**< Maximum number of iterations for the algorithm. */
    unsigned int currentNumIter;    /**< Current number of Picard algorithm iterations. */
};

/**
 * \struct TimeIncompressibleParams
 * \brief Time integration parameters.
 */
struct TimeIncompressibleParams
{
    bool adaptDT;               /**< Should the time step be changed during the computation? */
    double currentDT;           /**< Current time step. */
    double coeffDTincrease;     /**< Coefficient to multiply the current time step with if we increase it. */
    double coeffDTdecrease;     /**< Coefficient to divide the current time step with if we decrease it. */
    double maxDT;               /**< Maximu allowed time step. */
    double simuTime;            /**< Physical time which we want to simulate. */
    double simuDTToWrite;       /**< For which time interval should the program write data. */
    double nextWriteTrigger;    /**< Next physical time at which the program will write data. */
    double currentTime;
    unsigned int currentStep;
};

/**
 * \struct SolverIncompressibleParams
 * \brief Incompressible solver parameters.
 */
struct SolverIncompressibleParams
{
    double hchar;           /**< Characteristic size of an element. */
    double gravity;         /**< Acceleration of the gravity (g > 0). */
    FluidIncompressibleParams fluid;      /**< Fluid parameters. */
    PicardParams picard;    /**< Picard algorithm parameters. */
    TimeIncompressibleParams time;        /**< Time integration parameters. */
    std::array<bool, 5> whatToWrite {false}; /**< Which data will be written (u, v, p, ke, velocity). */
    std::string writeAs;
    std::string resultsName;    /**< File name in which the results will be written. */
    std::vector<double> initialCondition;    /**< Initial condition on u, v and p (mainly to have inlet) **/
#if defined(_OPENMP)
    unsigned int nOMPThreads;
#endif
};

/**
 * \class SolverIncompressible
 * \brief Represents a solver for an incompressible Newtonian fluid.
 */
class SolverIncompressible
{
    public:
        SolverIncompressible(const nlohmann::json& j, std::string mshName, std::string resultsName);
        ~SolverIncompressible();

        /**
         * \brief Display the parameters in SolverIncompressibleParams structure.
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
        SolverIncompressibleParams m_p;     /**< Solver parameters. */
        bool m_verboseOutput;               /**< Should the output be verbose? */

        Mesh m_mesh;                       /**< The mesh the solver is using. */

        std::vector<Eigen::MatrixXd> m_N;     /**< The shape functions matrices for eah gauss point (wihout *detJ). */
        Eigen::MatrixXd m_sumNTN;
        Eigen::VectorXd m_m;
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
};

#include "SolverIncompressible.inl"

#endif // SOLVERINCOMPRESSIBLE_HPP_INCLUDED
