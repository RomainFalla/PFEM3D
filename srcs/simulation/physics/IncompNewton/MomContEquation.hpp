#pragma once
#ifndef MOMCONTEQINCOMPNEWTON_HPP_INCLUDED
#define MOMCONTEQINCOMPNEWTON_HPP_INCLUDED

#include "../../Equation.hpp"
#include "../../nonLinearAlgo/PicardAlgo.hpp"

class Problem;
class Mesh;

class SIMULATION_API MomContEqIncompNewton : public Equation
{
    public:
        MomContEqIncompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                              std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                              const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex);
        MomContEqIncompNewton(const MomContEqIncompNewton& equation)             = delete;
        MomContEqIncompNewton& operator=(const MomContEqIncompNewton& equation)  = delete;
        MomContEqIncompNewton(MomContEqIncompNewton&& equation)                  = delete;
        MomContEqIncompNewton& operator=(MomContEqIncompNewton&& equation)       = delete;
        ~MomContEqIncompNewton() override;

        void displayParams() const override;

        bool solve() override;

    private:
        std::unique_ptr<PicardAlgo> m_pPicardAlgo;

        double m_rho;
        double m_mu;
        double m_gamma;
        double m_alpha;
        double m_Tr;

        Eigen::VectorXd m_bodyForce;

        void m_buildAb(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::VectorXd& qPrev);
        void m_applyBC(Eigen::VectorXd& b, const Eigen::VectorXd& qPrev);
        void m_executeTask(const Eigen::VectorXd& qPrev);

        double m_computeTauPSPG(const Element& element) const;
};

#endif // MOMCONTEQINCOMPNEWTON_HPP_INCLUDED
