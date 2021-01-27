#pragma once
#ifndef SOLVERWCOMPNEWTON_HPP_INCLUDED
#define SOLVERWCOMPNEWTON_HPP_INCLUDED

#include "../../Solver.hpp"

class SIMULATION_API SolverWCompNewton: public Solver
{
    public:
        SolverWCompNewton(Problem* pProblem, Mesh* pMesh, std::vector<SolTable> problemParams);
        SolverWCompNewton(const SolverWCompNewton& solver)             = delete;
        SolverWCompNewton& operator=(const SolverWCompNewton& solver)  = delete;
        SolverWCompNewton(SolverWCompNewton&& solver)                  = delete;
        SolverWCompNewton& operator=(SolverWCompNewton&& solver)       = delete;
        ~SolverWCompNewton() override;

        void displayParams() const override;

        bool solveOneTimeStep() override;
        void computeNextDT() override;

    protected:
        bool m_adaptDT;
		double m_securityCoeff;
		double m_maxDT;
		double m_initialDT;

		double m_nextTimeToRemesh;

		std::function<bool()> m_solveFunc;
		bool m_solveWCompNewtonNoT();
		bool m_solveBoussinesqWC();
};

#endif // SOLVERWCOMPNEWTON_HPP_INCLUDED
