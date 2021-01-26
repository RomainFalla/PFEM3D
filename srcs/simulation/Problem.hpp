#pragma once
#ifndef PROBLEM_HPP_INCLUDED
#define PROBLEM_HPP_INCLUDED

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "utility/SolTable.hpp"

#include "../mesh/Mesh.hpp"

#include "simulation_defines.h"

class Solver;
class Extractor;

/**
 * \class Problem
 * \brief Represents a problem which will be solved by a solver.
 *
 * The Problem class is responsible for holding general informations about the current
 * problem being simulated: the end simulation time, the current time and step, the
 * number of threads used or which and how data can be written.
 */

class SIMULATION_API Problem
{
    public:
        Problem(const std::string& luaFilePath);
        Problem(const Problem& problem)             = delete;
        Problem& operator=(const Problem& problem)  = delete;
        Problem(Problem&& problem)                  = delete;
        Problem& operator=(Problem&& problem)       = delete;
        virtual ~Problem();

        /// \brief Display general parameters of the problem.
        virtual void displayParams() const;

        /// \return The id of the problem (child class have to set m_id).
        std::string getID() const noexcept;

        /// \return A vector containing the name of the data which can be extracted for this problem.
        virtual std::vector<std::string> getWrittableDataName() const;

        /// \param name The name of the data to write.
        /// \param nodeIndex From which node should the data be extracted.
        /// \return A vector containing the value of the data.
        virtual std::vector<double> getWrittableData(const std::string& name, std::size_t nodeIndex) const;

        /// \return A vector containing the name of the global data which can be extracted for this problem.
        virtual std::vector<std::string> getGlobalWrittableDataName() const;

        /// \param name The name of the global data to write.
        /// \return The value of the global data.
        virtual double getGlobalWrittableData(const std::string& name) const;

        /// \return A vector containing the name of the mesh data which can be extracted for this problem.
        std::vector<std::string> getMeshWrittableDataName() const;

        /// \param name The name of the mesh data to write.
        /// \return The value of the mesh data.
        std::vector<double> getMeshWrittableData(const std::string& name, std::size_t nodeIndex) const;

        inline const Mesh& getMesh() const noexcept;
        inline const Solver& getSolver() const noexcept;
        inline double getCurrentSimTime() const noexcept;
        inline double getCurrentSimStep() const noexcept;
        inline unsigned int getThreadCount() const noexcept;
        inline bool isOutputVerbose() const noexcept;

        /// \brief This function simulate the problem from t = 0 to t = t_max
        void simulate();

        /// \param timeStep The time step
        /// \brief Update the internal current time with the time step
        void updateTime(double timeStep);

    protected:
        std::string m_id;               /**<  The id of the problem (should be set by child class). */
        double m_time;                  /**<  Current time of the simulation. */
        double m_maxTime;               /**<  Maximum time of the simulation. */
        std::size_t m_step;             /**<  Current step of the simulation. */
        unsigned int m_statesNumber;    /**<  Number of states used by the problem. */

        bool m_verboseOutput;           /**<  Should the console output be verbose ?. */

        unsigned int m_nThreads;        /**<  The number of OpenMP threads used. */
        std::vector<sol::state> m_solState;     /**<  sol::state's for the lua parameters file (1 per thread). */
        std::vector<SolTable> m_problemParams;  /**<  sol::table wrapper for the problem parameters (1 per thread). */

        std::unique_ptr<Mesh> m_pMesh;          /**<  Smart pointer to the mesh. */
        std::unique_ptr<Solver> m_pSolver;      /**<  Smart pointer to the solver. */
        std::vector<std::unique_ptr<Extractor>> m_pExtractors;  /**<  Smart pointer to the extractirs. */

        /// \brief This function parse the lua parameters file and set the requires extractors in m_pExtractors.
        /// Should be called in the constructor of every child class.
        void addExtractors();

        /// \brief This function parse the lua parameters file and set the initial condition.
        /// Should be called in the constructor of every child class.
        void setInitialCondition();
};

#include "Problem.inl"

#endif // PROBLEM_HPP_INCLUDED
