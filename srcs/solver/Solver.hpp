#pragma once
#ifndef SOLVER_HPP_INCLUDED
#define SOLVER_HPP_INCLUDED

#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifndef SOL_ALL_SAFETIES_ON
    #define SOL_ALL_SAFETIES_ON 1
#endif
#include <sol3/sol.hpp>

#if defined(_OPENMP)
    #include <omp.h>
#endif // defined

#include "../mesh/Mesh.hpp"
#include "extractors/Extractors.hpp"

#include "Solver_export.h"

enum class SOLVER_TYPE
{
    Undefined,
    Incompressible_PSPG,
    WeaklyCompressible
};

struct SolverCreateInfo
{
    MeshCreateInfo meshInfos = {};
    double gravity = 0;
    bool strongPAtFS = false;
    bool adaptDT = true;
    double initialDT = 0;
    double endTime = 0;
    double maxDT = 0;
    std::string IBCfile = {};
};

/**
 * \class Solver
 * \brief Represents a solver.
 */
class SOLVER_API Solver
{
    public:
        Solver()                                = delete;
        Solver(const SolverCreateInfo& solverInfos);
        Solver(const Solver& solver)            = delete;
        Solver& operator=(const Solver& solver) = delete;
        Solver(Solver&& solver)                 = delete;
        Solver& operator=(Solver&& solver)      = delete;
        virtual ~Solver() noexcept              = default;

        void addPointExtractor(const std::string& outFileName, double timeBetweenWriting,
                               unsigned short stateToWrite, const std::vector<std::vector<double>>& points);

        void addMinMaxExtractor(const std::string& outFileName, double timeBetweenWriting,
                                unsigned short coordinate, const std::string& minMax);

        void addMassExtractor(const std::string& outFileName, double timeBetweenWriting);

        void addGMSHExtractor(const std::string& outFileName, double timeBetweenWriting,
                              const std::vector<std::string>& whatToWrite, std::string writeAs);

        /**
         * \return Return the current time increment of the simulation in seconds.
         */
        inline double getCurrentDT() const noexcept;

        /**
         * \return Return the current step of the simulation.
         */
        inline unsigned int getCurrentStep() const noexcept;

        /**
         * \return Return the current physical time of the simulation in seconds.
         */
        inline double getCurrentTime() const noexcept;

        /**
         * \param elm The element index.
         * \param state The considered state.
         * \return Return a vector containing the state of a whole element.
         */
        inline Eigen::VectorXd getElementState(IndexType elm, unsigned short state) const noexcept;

        /**
         * \return Return the maximum physical time at which the simulation will stop in seconds.
         */
        inline double getEndTime() const noexcept;

        /**
         * \return Return the maximum time increment authorized in seconds.
         */
        inline double getMaxDT() const noexcept;

        /**
         * \return Return a constant reference to the mesh used by the solver.
         */
        inline const Mesh& getMesh() const noexcept;

        /**
         * \return A Eigen vector containing the states of the node.
         * \param beginState The first state which will be contained in q.
         * \param endState The last state which will be contained in q.
         */
        inline Eigen::VectorXd getQFromNodesStates(unsigned short beginState, unsigned short endState) const noexcept;

        /**
        * \return What solver it is.
        */
        inline SOLVER_TYPE getSolverType() const noexcept;

        /**
        * \return Return the number of state the solver is using.
        */
        inline unsigned short getStatesNumber() const noexcept;

        /**
         * \return Return if the current time step can be modified or not.
         */
        inline bool isDTAdaptable() const noexcept;

        /**
         * \brief Update the time interval for the next time step (only if authorized).
         * \param dt The new time interval.
         */
        inline void setCurrentDT(double dt);

        /**
         * \brief Set the states of the nodes from a Eigen vector.
         * \param q a Eigen vector containing the value of the new states. Its size is (endState-beginState+1).
         * \param beginState The first state contained in q.
         * \param endState The last state contained in q.
         */
        inline void setNodesStatesfromQ(const Eigen::VectorXd& q, unsigned short beginState, unsigned short endState) noexcept;

    protected:
        SOLVER_TYPE m_solverType;

        sol::state m_lua;

        unsigned int m_numOMPThreads;   /**< Number of OpenMP threads used. */

        unsigned short m_statesNumber;

        double m_gravity;               /**< Acceleration of the gravity (g > 0). */
        bool m_strongPAtFS;

        //Time Parameters
        bool m_adaptDT;                 /**< Should the time step be changed during the computation? */
        double m_currentDT;             /**< Current time step. */
        unsigned int m_currentStep;
        double m_currentTime;
        double m_endTime;                /**< How many real seconds do we compute ? */
        double m_maxDT;                 /**< Maximum allowed time step. */

        Mesh m_mesh;                       /**< The mesh the solver is using. */

        std::vector<Eigen::MatrixXd> m_N;     /**< The shape functions matrices for eah gauss point (wihout *detJ). */
        Eigen::VectorXd m_m;
        Eigen::VectorXd m_bodyForces;

        std::vector<std::unique_ptr<Extractor>> m_pExtractor;
        bool m_hasGMSHExtractor;

        /**
         * \param elementIndex The index of the element in the element list.
         * \return The gradient shape function matrix for the element in the format:
         *         [dN1dx dN2dx dN3dx 0 0 0; 0 0 0 dN1dy dN2dy dN3dy; dN1dx dN2dx dN3dx dN1dy dN2dy dN3dy]
         */
        inline Eigen::MatrixXd getB(IndexType elementIndex) const noexcept;

        /**
         * \return The gradient shape function matrix for each gauss point in the format:
         *         [N1 N2 N3 0 0 0; 0 0 0 N1 N2 N3]
         */
        inline std::vector<Eigen::MatrixXd> getN() const noexcept;

        /**
         * \brief Set the initial condition on u, v, p for the initial cloud of nodes.
         */
        void setInitialCondition();
};

#include "Solver.inl"

#endif // SOLVER_HPP_INCLUDED
