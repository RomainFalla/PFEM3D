#ifndef SOLVER_HPP_INCLUDED
#define SOLVER_HPP_INCLUDED

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../params/Params.hpp"
#include "../mesh/Mesh.hpp"

enum ProblemType
{
    INCOMPRESSIBLE_PSPG
};

class Solver
{
    public:
        Solver(const Params& params, Mesh& mesh, ProblemType ProblemType);
        ~Solver();

        void solveProblem();

    private:
        std::vector<bool> _getIndices() const;
        void _buildPicardSystem();
        void _computeTauPSPG();
        bool _solveSystem();

        const Params& m_params;
        Mesh& m_mesh;

        ProblemType m_problemType;

        Eigen::VectorXd m_qprev;
        Eigen::VectorXd m_q;
        Eigen::SparseMatrix<double> m_A;
        Eigen::VectorXd m_b;

        Eigen::SparseMatrix<double> m_M;
        Eigen::SparseMatrix<double> m_K;
        Eigen::SparseMatrix<double> m_D;
        Eigen::SparseMatrix<double> m_C;
        Eigen::SparseMatrix<double> m_L;
        Eigen::VectorXd m_F;
        Eigen::VectorXd m_H;

        std::vector<double> m_tauPSPG;

        double m_currentDT;

        double m_numIterSolverMsh;

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solver;
};

#endif // SOLVER_HPP_INCLUDED
