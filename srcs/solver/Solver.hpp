#ifndef SOLVER_HPP_INCLUDED
#define SOLVER_HPP_INCLUDED

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../params/Params.hpp"
#include "../mesh/Mesh.hpp"

/**
 * \struct FluidParams
 * \brief Incompressible fluid parameters.
 */
struct FluidParams
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
 * \struct TimeParams
 * \brief Time integration parameters.
 */
struct TimeParams
{
    bool adaptDT;               /**< Should the time step be changed during the computation? */
    double currentDT;           /**< Current time step. */
    double coeffDTincrease;     /**< Coefficient to multiply the current time step with if we increase it. */
    double coeffDTdecrease;     /**< Coefficient to divide the current time step with if we decrease it. */
    double maxDT;               /**< Maximu allowed time step. */
    double simuTime;            /**< Physical time which we want to simulate. */
    double simuDTToWrite;       /**< For which time interval should the program write data. */
    double nextWriteTrigger;    /**< Next physical time at which the program will write data. */
};

/**
 * \struct SolverIncompressibleParams
 * \brief Incompressible solver parameters.
 */
struct SolverIncompressibleParams
{
    double hchar;           /**< Characteristic size of an element. */
    double gravity;         /**< Acceleration of the gravity (g > 0). */
    FluidParams fluid;      /**< Fluid parameters. */
    PicardParams picard;    /**< Picard algorithm parameters. */
    TimeParams time;        /**< Time integration parameters. */
    std::array<bool, 5> whatToWrite {false}; /**< Which data will be written (u, v, p, ke, velocity). */
    std::string resultsName;    /**< File name in which the results will be written. */
    std::vector<double> initialCondition;    /**< Initial condition on u, v and p (mainly to have inlet) **/
};

/**
 * \class Solver
 * \brief Represents a solver for an incompressible Newtonian fluid.
 */
class Solver
{
    public:
        Solver(const Params& params, Mesh& mesh, std::string resultsName);
        ~Solver();

        /**
         * \brief Solve the problem for a certain set of parameters.
         */
        void solveProblem();

    private:
        /**
         * \brief Return a vector of boolean to know which line of the A matrix should
         * transformed into ... 0 1 0 ... and which element of the b vector set to qprev.
         * \return Vector of boolean in which an element is true if the node is a
         *         boundary node (u and v should be conserved) or the node is a free
         *         node (u, v, p) should be conservec.
         */
        std::vector<bool> _getIndices() const;

        /**
         * \brief Set the initial condition on u, v, p for the initial cloud of nodes.
         */
        void _setInitialCondition();

        /**
         * \brief Build the matrix A and the vector b of the Picard Algorithm.
         */
        void _buildPicardSystem();

        /**
         * \brief Build the matrix A and the vector b of the Picard Algorithm.
         */
        void _computeTauPSPG();

        /**
         * \brief Solve the Picard algorithm for one time step.
         * \return true if the algorithm converged, false otherwise.
         */
        bool _solveSystem();

        /**
         * \brief Write data to a file.
         * \param time Current physical time
         * \param step Current step.
         */
        void _write(double time, unsigned int step);
        //void _writeFinalize();

        Mesh& m_mesh;                       /**< The mesh the solver is using. */

        Eigen::VectorXd m_qprev;            /**< The precedent solution */
        Eigen::SparseMatrix<double> m_A;    /**< The matrix A representing the problem:
                                                 [M/dt+K -D^T; C/dt-D L]. */
        Eigen::VectorXd m_b;                /**< The vector b representing the problem:
                                                 [M/dt*qprev + F; H]. */

        Eigen::SparseMatrix<double> m_M;     /**< The mass matrix. */
        Eigen::SparseMatrix<double> m_K;     /**< The viscosity matrix. */
        Eigen::SparseMatrix<double> m_D;     /**< The pressure matrix. */
        Eigen::VectorXd m_F;                 /**< The volume force vector. */

        std::vector<double> m_tauPSPG;       /**< tau_PSPG parameters for each element. */

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solverLU; /**< Eigen SparseLU solver. */

        SolverIncompressibleParams m_p;     /**< Solver parameters. */

        bool m_verboseOutput;               /**< Should the output be verbose? */
};

#endif // SOLVER_HPP_INCLUDED
