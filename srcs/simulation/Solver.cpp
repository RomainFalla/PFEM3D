#include "Solver.hpp"

#include "../mesh/Mesh.hpp"
#include "Problem.hpp"
#include "Equation.hpp"

Solver::Solver(Problem* pProblem, Mesh* pMesh, std::vector<SolTable> m_problemParams):
m_timeStep(0),
m_pMesh(pMesh),
m_pProblem(pProblem)
{
    m_solverParams.resize(m_pProblem->getThreadCount());
    for(std::size_t i = 0 ; i < m_solverParams.size() ; ++i)
        m_solverParams[i] = SolTable("Solver", m_problemParams[i]);

    m_id = m_solverParams[0].checkAndGet<std::string>("id");
}

Solver::~Solver()
{

}

void Solver::displayParams() const
{
    throw std::runtime_error("Unimplemented function by the child class -> Solver::displayParams()");
}

std::string Solver::getID() const noexcept
{
    return m_id;
}

bool Solver::solveOneTimeStep()
{
    throw std::runtime_error("Unimplemented function by the child class -> Solver::solveOneTimeStep()");
}

void Solver::computeNextDT()
{
    throw std::runtime_error("Unimplemented function by the child class -> Solver::computeNextDT()");
}

bool Solver::checkBC(SolTable bcParam, unsigned int n, const Node& node, std::string bcString, unsigned int expectedBCSize)
{
    bool res = bcParam.checkCallNoThrow(m_pMesh->getNodeType(n) + bcString,
                                        node.getPosition(),
                                        m_pMesh->getBoundNodeInitPos(n),
                                        node.getStates(), 0);

    if(res)
    {
        std::vector<double> result = bcParam.call<std::vector<double>>(
            m_pMesh->getNodeType(n) + bcString,
            node.getPosition(),
            m_pMesh->getBoundNodeInitPos(n),
            node.getStates(),
            0);

        if(result.size() != expectedBCSize)
        {
            throw std::runtime_error(" the boundary condition " + m_pMesh->getNodeType(n) + bcString +
                                     " does not return the right number of states: " +
                                     std::to_string(result.size()) + " vs " +  std::to_string(expectedBCSize));
        }
    }

    return res;
}

