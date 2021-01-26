#pragma once
#ifndef EQUATION_HPP_INCLUDED
#define EQUATION_HPP_INCLUDED

#include <vector>
#include <Eigen/Dense>

#include "utility/SolTable.hpp"
#include "matricesBuilder/MatricesBuilder.hpp"

#include "simulation_defines.h"

class Problem;
class Solver;
class Mesh;
class Node;

/**
 * \class Equation
 * \brief Represents one equation of a problem.
 *
 * This class is responsible of solving one equation and update the corresponding
 * mesh states accordingly.
 */
class SIMULATION_API Equation
{
    public:
        Equation(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                 std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                 const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex,
                 const std::string& id);
        Equation(const Equation& equation)             = delete;
        Equation& operator=(const Solver& equation)    = delete;
        Equation(Equation&& equation)                  = delete;
        Equation& operator=(Equation&& equation)       = delete;
        virtual ~Equation();

        /// \brief Display general parameters of the equation.
        virtual void displayParams() const;

        /// \param thread the OpenMP thread which might be calling this function.
        /// \return The underlying boundary conditions parameters.
        SolTable getBCParam(unsigned int thread) const;

        /// \return The id of the equation (which will be used to find specific equation parameters in the lua file).
        std::string getID() const noexcept;

        /// \return Does this equation requires the computation of normals and curvatures (child class should set m_needNormalCurv) ?
        bool isNormalCurvNeeded() const noexcept;

        /// \param he The element characteristic size.
        /// \param node The node on which the speed is calculated.
        /// \return An equivalent speed if the problem is solved using an explicit method.
        virtual double getSpeedEquiv(double he, const Node& node);

        /// \brief A function to precompute data for the equation of needed.
        virtual void preCompute();

        /// \brief Solve the equation and update the mesh states accordingly.
        virtual bool solve();


    protected:
        bool m_needNormalCurv;  /**< Does this equation require the computation of the normals and curvatures ? */
        std::string m_id;       /**<  The id of the equation. */
        std::vector<SolTable> m_materialParams; /**<  sol::table wrapper for the material parameters (1 per thread). */
        std::vector<SolTable> m_equationParams; /**<  sol::table wrapper for the equation parameters (1 per thread). */
        std::vector<SolTable> m_bcParams;       /**<  sol::table wrapper for the boundary condition parameters (1 per thread). */

        std::vector<unsigned short> m_bcFlags;  /**< Flags that can be used by the equation to know which boundary condition is set for a node. */
        std::vector<unsigned int> m_statesIndex;    /**< Index of the states set by this equation in the mesh */

        Problem* m_pProblem;    /**< Pointer to the underlying problem. */
        Solver* m_pSolver;      /**< Pointer to the underlying solver. */
        Mesh* m_pMesh;          /**< Pointer to the mesh used. */

        std::unique_ptr<MatrixBuilder> m_pMatBuilder; /**< Class responsible of building the required matrices. */
};

#endif // EQUATION_HPP_INCLUDED
