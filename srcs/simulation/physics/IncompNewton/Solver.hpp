#pragma once
#ifndef SOLVERINCOMPNEWTON_HPP_INCLUDED
#define SOLVERINCOMPNEWTON_HPP_INCLUDED

#include "../../Solver.hpp"

class SIMULATION_API SolverIncompNewton: public Solver
{
    public:
        SolverIncompNewton(Problem* pProblem, Mesh* pMesh, std::vector<SolTable> problemParams);
        SolverIncompNewton(const SolverIncompNewton& solver)             = delete;
        SolverIncompNewton& operator=(const SolverIncompNewton& solver)  = delete;
        SolverIncompNewton(SolverIncompNewton&& solver)                  = delete;
        SolverIncompNewton& operator=(SolverIncompNewton&& solver)       = delete;
        ~SolverIncompNewton() override;

        void displayParams() const override;

        bool solveOneTimeStep() override;
        void computeNextDT() override;

    protected:
        bool m_adaptDT;
		double m_coeffDTincrease;
		double m_coeffDTDecrease;
		double m_maxDT;
		double m_initialDT;

		bool m_solveHeatFirst;

		std::function<bool()> m_solveFunc;
		bool m_solveIncompNewtonNoT();
		bool m_solveBoussinesq();
		bool m_solveConduction();

};

#endif // SOLVERINCOMPNEWTON_HPP_INCLUDED
