#pragma once
#ifndef HEATEQINCOMPNEWTON_HPP_INCLUDED
#define HEATEQINCOMPNEWTON_HPP_INCLUDED

#include "../../Equation.hpp"
#include "../../nonLinearAlgo/PicardAlgo.hpp"

class Problem;
class Mesh;

class SIMULATION_API HeatEqIncompNewton : public Equation
{
    public:
        HeatEqIncompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                          std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                          const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex);
        HeatEqIncompNewton(const HeatEqIncompNewton& equation)             = delete;
        HeatEqIncompNewton& operator=(const HeatEqIncompNewton& equation)  = delete;
        HeatEqIncompNewton(HeatEqIncompNewton&& equation)                  = delete;
        HeatEqIncompNewton& operator=(HeatEqIncompNewton&& equation)       = delete;
        ~HeatEqIncompNewton() override;

        void displayParams() const override;

        bool solve() override;

    private:
        std::unique_ptr<PicardAlgo> m_pPicardAlgo;

        double m_k;
        double m_cv;
        double m_rho;

        void m_buildAb(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::VectorXd& qPrev);
        void m_applyBC(Eigen::VectorXd& b, const Eigen::VectorXd& qPrev);
        void m_executeTask(const Eigen::VectorXd& qPrev);

        double m_computeTauPSPG(const Element& element) const;
};

#endif // HEATEQINCOMPNEWTON_HPP_INCLUDED
