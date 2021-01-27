#pragma once
#ifndef PROBWCOMPNEWTON_HPP_INCLUDED
#define PROBWCOMPNEWTON_HPP_INCLUDED

#include "../../Problem.hpp"

/** \class ProbWCompNewton
 * Represent a weakly-compressible Newtonian flow with Tait-Murnaghan state equation  without temperature:
 * \f{eqnarray*}{
 *   \begin{cases}
 *   \dfrac{d}{dt}\rho + \rho\vec{\nabla}\cdot\vec{u} = 0 & \forall \vec{x}\in V\\
 *   \dfrac{d}{dt}\vec{v} = \vec{\nabla}\cdot\vec{\vec{\sigma}} + \rho\vec{b} & \forall \vec{x}\in V\\
 *   \vec{\vec{\sigma}} = -p\vec{\vec{I}} + 2\mu\textrm{dev}\vec{\vec{D}} & \forall \vec{x}\in V\\
 *   p = \dfrac{K_0}{K_0^\prime}\left[\left(\dfrac{\rho}{\rho^\star}\right)^{K_0^\prime}-1\right] & \forall \vec{x}\in V\\
 *   \vec{v}(\vec{x}, t=0) = \vec{v}_0,\quad p(\vec{x}, t=0) = p_0 ,\quad \rho(\vec{x}, t=0) = \rho_0 & \forall \vec{x}\in V_0\\
 *   \vec{v} = \vec{v}_b & \forall \vec{x}\in S_D\\
 *   \vec{\vec{\sigma}}\cdot\hat{n} = \gamma\kappa\hat{n} & \forall \vec{x}\in S_N
 *   \end{cases}
 * \f}
 */
class SIMULATION_API ProbWCompNewton : public Problem
{
    public:
        ProbWCompNewton(std::string luaFilePath);
        ProbWCompNewton(const ProbWCompNewton& problem)             = delete;
        ProbWCompNewton& operator=(const ProbWCompNewton& problem)  = delete;
        ProbWCompNewton(ProbWCompNewton&& problem)                  = delete;
        ProbWCompNewton& operator=(ProbWCompNewton&& problem)       = delete;
        ~ProbWCompNewton() override;

        void displayParams() const override;

        std::vector<std::string> getWrittableDataName() const override;
        std::vector<double> getWrittableData(const std::string& name, std::size_t nodeIndex) const override;
        std::vector<std::string> getGlobalWrittableDataName() const override;
        double getGlobalWrittableData(const std::string& name) const override;

    private:

};

#endif // PROBWCOMPNEWTON_HPP_INCLUDED
