#pragma once
#ifndef MOMEQWCOMPNEWTON_HPP_INCLUDED
#define MOMEQWCOMPNEWTON_HPP_INCLUDED

#include "../../Equation.hpp"

class Problem;
class Mesh;

class SIMULATION_API MomEqWCompNewton : public Equation
{
    public:
        MomEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                          std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                          const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex);
        MomEqWCompNewton(const MomEqWCompNewton& equation)             = delete;
        MomEqWCompNewton& operator=(const MomEqWCompNewton& equation)  = delete;
        MomEqWCompNewton(MomEqWCompNewton&& equation)                  = delete;
        MomEqWCompNewton& operator=(MomEqWCompNewton&& equation)       = delete;
        ~MomEqWCompNewton() override;

        void displayParams() const override;

        double getSpeedEquiv(double he, const Node& node) override;
        bool solve() override;
        void setQVhalf(Eigen::VectorXd qV1half);

    private:
        double m_mu;
        double m_gamma;
        double m_alpha;
        double m_Tr;

        Eigen::VectorXd m_bodyForce;

        void m_buildSystem(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& A, Eigen::VectorXd& b);
        void m_applyBC(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& A, Eigen::VectorXd& b);
};

#endif // MOMEQWCOMPNEWTON_HPP_INCLUDED
