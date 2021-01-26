#include "Equation.hpp"
#include "Problem.hpp"

Equation::Equation(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                 std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                 const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex,
                 const std::string& id):
m_id(id),
m_materialParams(materialParams),
m_bcFlags(bcFlags),
m_statesIndex(statesIndex),
m_pProblem(pProblem),
m_pSolver(pSolver),
m_pMesh(pMesh)
{
    unsigned int nGPHD = 0;
    unsigned int nGPLD = 0;
    if(m_pMesh->getDim() == 2)
    {
        nGPHD = 3;
        nGPLD = 3;
    }
    else
    {
        nGPHD = 4;
        nGPLD = 3;
    }
    m_pMatBuilder = std::make_unique<MatrixBuilder>(*pMesh, nGPHD, nGPLD);

    m_equationParams.resize(m_pProblem->getThreadCount());
    m_bcParams.resize(m_pProblem->getThreadCount());
    for(std::size_t i = 0 ; i < solverParams.size() ; ++i)
    {
        m_equationParams[i] = SolTable(getID(), solverParams[i]);
        m_bcParams[i] = SolTable("BC", m_equationParams[i]);
    }
}

Equation::~Equation()
{

}

void Equation::displayParams() const
{
    throw std::runtime_error("Unimplemented function by the child class -> Equation::displayParams()");
}

SolTable Equation::getBCParam(unsigned int thread) const
{
    return m_bcParams[thread];
}

std::string Equation::getID() const noexcept
{
    return m_id;
}

bool Equation::isNormalCurvNeeded() const noexcept
{
    return m_needNormalCurv;
}

double Equation::getSpeedEquiv(double /** he **/, const Node& /** node **/)
{
    throw std::runtime_error("Unimplemented function by the child class -> Equation::getSpeedEquiv()!");
}

void Equation::preCompute()
{
    throw std::runtime_error("Unimplemented function by the child class -> Equation::preCompute()!");
}

bool Equation::solve()
{
    throw std::runtime_error("Unimplemented function by the child class -> Equation::solve()!");
}
