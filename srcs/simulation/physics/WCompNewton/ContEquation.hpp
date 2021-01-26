#pragma once
#ifndef CONTEQWCOMPNEWTON_HPP_INCLUDED
#define CONTEQWCOMPNEWTON_HPP_INCLUDED

#include "../../Equation.hpp"

class Problem;
class Mesh;

class SIMULATION_API ContEqWCompNewton : public Equation
{
    public:
        ContEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                          std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                          const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex);
        ContEqWCompNewton(const ContEqWCompNewton& equation)             = delete;
        ContEqWCompNewton& operator=(const ContEqWCompNewton& equation)  = delete;
        ContEqWCompNewton(ContEqWCompNewton&& equation)                  = delete;
        ContEqWCompNewton& operator=(ContEqWCompNewton&& equation)       = delete;
        ~ContEqWCompNewton() override;

        void displayParams() const override;

        double getSpeedEquiv(double he, const Node& node) override;
        bool solve() override;
        void preCompute() override;

    private:
        double m_K0;
        double m_K0p;
        double m_rhoStar;
        bool m_strongContinuity;
        double m_keStab;

        Eigen::VectorXd m_F0;

        void m_buildF0(Eigen::VectorXd& F0);
        void m_buildSystem(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F);
        void m_applyBC(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F);
        Eigen::VectorXd m_getPFromRhoTaitMurnagham(const Eigen::VectorXd& rho);
};

#endif // CONTEQWCOMPNEWTON_HPP_INCLUDED
