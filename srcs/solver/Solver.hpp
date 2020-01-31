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
        void _buildMatrices(std::vector<bool> matricesToBuild, std::vector<bool> vectorsToBuild);
        void _computeSideVariables();

        const Params& params;
        Mesh& mesh;

        ProblemType problemType;

        std::vector<Eigen::VectorXd> q;
        std::vector<Eigen::SparseMatrix<double>> matrices;
        std::vector<Eigen::VectorXd> vectors;

        std::vector<std::vector<double>> sideVariables; //HelloTauPSPG

        double currentTimeStep;

};
