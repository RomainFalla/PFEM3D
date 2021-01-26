#pragma once
#ifndef PROBINCOMPNEWTON_HPP_INCLUDED
#define PROBINCOMPNEWTON_HPP_INCLUDED

#include "../../Problem.hpp"

/** \class ProbIncompNewton
 * Represent an incompressible Newtonian flow without temperature:
 * \f{eqnarray*}{
 *   \begin{cases}
 *   \vec{\nabla}\cdot\vec{u} = 0 & \forall \vec{x}\in V\\
 *   \dfrac{d}{dt}\vec{v} = \vec{\nabla}\cdot\vec{\vec{\sigma}} + \rho\vec{b} & \forall \vec{x}\in V\\
 *   \vec{\vec{\sigma}} = -p\vec{\vec{I}} + 2\mu\vec{\vec{D}} & \forall \vec{x}\in V\\
 *   \vec{v}(\vec{x}, t=0) = \vec{v}_0,\quad p(\vec{x}, t=0) = p_0  & \forall \vec{x}\in V_0\\
 *   \vec{v} = \vec{v}_b & \forall \vec{x}\in S_D\\
 *   \vec{\vec{\sigma}}\cdot\hat{n} = \gamma\kappa\hat{n} & \forall \vec{x}\in S_N
 *   \end{cases}
 * \f}
 */
class SIMULATION_API ProbIncompNewton : public Problem
{
    public:
        ProbIncompNewton(std::string luaFilePath);
        ProbIncompNewton(const ProbIncompNewton& problem)             = delete;
        ProbIncompNewton& operator=(const ProbIncompNewton& problem)  = delete;
        ProbIncompNewton(ProbIncompNewton&& problem)                  = delete;
        ProbIncompNewton& operator=(ProbIncompNewton&& problem)       = delete;
        ~ProbIncompNewton() override;

        void displayParams() const override;

        std::vector<std::string> getWrittableDataName() const override;
        std::vector<double> getWrittableData(const std::string& name, std::size_t nodeIndex) const override;
        std::vector<std::string> getGlobalWrittableDataName() const override;
        double getGlobalWrittableData(const std::string& name) const override;

    private:

};

#endif // PROBINCOMPNEWTON_HPP_INCLUDED
