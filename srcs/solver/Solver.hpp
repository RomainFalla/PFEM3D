#ifndef SOLVER_HPP_INCLUDED
#define SOLVER_HPP_INCLUDED

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../params/Params.hpp"
#include "../mesh/Mesh.hpp"

struct FluidParams
{
    double rho;
    double mu;
};

struct PicardParams
{
    double relTol;
    unsigned int maxIter;
    unsigned int currentNumIter;
};

struct TimeParams
{
    bool adaptDT;
    double currentDT;
    double coeffDTincrease;
    double coeffDTdecrease;
    double maxDT;
    double simuTime;
    double simuDTToWrite;
    double nextWriteTrigger;
};

struct SolverIncompressibleParams
{
    double hchar;
    double gravity;
    FluidParams fluid;
    PicardParams picard;
    TimeParams time;
    std::array<bool, 5> whatToWrite {false}; //u v p ke velocity
    std::string resultsName;
};

class Solver
{
    public:
        Solver(const Params& params, Mesh& mesh, std::string resultsName);
        ~Solver();

        void solveProblem();

    private:
        std::vector<bool> _getIndices() const;
        void _buildPicardSystem();
        void _computeTauPSPG();
        bool _solveSystem();
        void _write(double time, unsigned int step);
        void _writeFinalize();

        Mesh& m_mesh;

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

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solver;

        SolverIncompressibleParams m_p;

        bool m_verboseOutput;
};

#endif // SOLVER_HPP_INCLUDED
