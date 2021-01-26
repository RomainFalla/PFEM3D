#pragma once
#ifndef HEATEQWCOMPNEWTON_HPP_INCLUDED
#define HEATEQWCOMPNEWTON_HPP_INCLUDED

#include "../../Equation.hpp"

class Problem;
class Mesh;

class SIMULATION_API HeatEqWCompNewton : public Equation
{
    public:
        HeatEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                          std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                          const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex);
        HeatEqWCompNewton(const HeatEqWCompNewton& equation)             = delete;
        HeatEqWCompNewton& operator=(const HeatEqWCompNewton& equation)  = delete;
        HeatEqWCompNewton(HeatEqWCompNewton&& equation)                  = delete;
        HeatEqWCompNewton& operator=(HeatEqWCompNewton&& equation)       = delete;
        ~HeatEqWCompNewton() override;

        void displayParams() const override;

        double getSpeedEquiv(double he, const Node& node) override;
        bool solve() override;

    private:
        double m_k;
        double m_cv;

        void m_buildSystem(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F);
        void m_applyBC(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F);
};

#endif // HEATEQWCOMPNEWTON_HPP_INCLUDED
