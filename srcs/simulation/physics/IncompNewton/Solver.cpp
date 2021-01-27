#include "Problem.hpp"
#include "Solver.hpp"
#include "MomContEquation.hpp"
#include "HeatEquation.hpp"


SolverIncompNewton::SolverIncompNewton(Problem* pProblem, Mesh* pMesh, std::vector<SolTable> problemParams):
Solver(pProblem, pMesh, problemParams)
{
    //Check if the asked problem and solver are supported
    if(m_pProblem->getID() != "IncompNewtonNoT" && m_pProblem->getID() != "Boussinesq" && m_pProblem->getID() != "Conduction")
        throw std::runtime_error("this solver cannot be used with problem whose id is " + m_pProblem->getID());

    if(m_id != "PSPG")
        throw std::runtime_error("this solver does not know id " + m_id);

    //Load material params for equations
    std::vector<SolTable> materialParams(problemParams.size());
    for(std::size_t i = 0 ; i < problemParams.size() ; ++i)
        materialParams[i] = SolTable("Material", problemParams[i]);

    //Load equations depending of the problem
    std::vector<unsigned short> bcFlags;
    std::vector<unsigned int> statesIndex;
    if(m_pProblem->getID() == "IncompNewtonNoT")
    {
        //Only the momentum-continuity equation
        m_pEquations.resize(1);

        bcFlags = {0};
        statesIndex = {0};
        m_pEquations[0] = std::make_unique<MomContEqIncompNewton>(
            m_pProblem, this, m_pMesh, m_solverParams, materialParams,
            bcFlags, statesIndex
        );

        //Set the right node flag if the boundary condition is present
        SolTable bcParam = m_pEquations[0]->getBCParam(0);
        for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
        {
            const Node& node = m_pMesh->getNode(n);
            if(node.isBound())
            {
                bool res = checkBC(bcParam, n, node, "V", m_pMesh->getDim());

                if(res)
                    m_pMesh->setNodeFlag(n, 0);
            }
        }

        m_solveFunc = std::bind(&SolverIncompNewton::m_solveIncompNewtonNoT, this);
    }
    else if(m_pProblem->getID() == "Boussinesq")
    {
        //The momentum-continuity equation and the heat equation
        m_pEquations.resize(2);

        bcFlags = {0};
        statesIndex = {0, static_cast<unsigned int>(m_pMesh->getDim()) + 1};
        m_pEquations[0] = std::make_unique<MomContEqIncompNewton>(
            m_pProblem, this, m_pMesh, m_solverParams, materialParams,
            bcFlags, statesIndex
        );

        bcFlags = {1, 2}; //Dirichlet and Neumann
        statesIndex = {static_cast<unsigned int>(m_pMesh->getDim()) + 1};
        m_pEquations[1] = std::make_unique<HeatEqIncompNewton>(
            m_pProblem, this, m_pMesh, m_solverParams, materialParams,
            bcFlags, statesIndex
        );

        //Set the right node flag if the boundary condition is present
        SolTable bcParamMomCont = m_pEquations[0]->getBCParam(0);
        SolTable bcParamHeat = m_pEquations[1]->getBCParam(0);
        for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
        {
            const Node& node = m_pMesh->getNode(n);
            if(node.isBound())
            {
                bool resV = checkBC(bcParamMomCont, n, node, "V", m_pMesh->getDim());
                bool resT = checkBC(bcParamHeat, n, node, "T", 1);
                bool resQ = checkBC(bcParamHeat, n, node, "Q", 1);

                if(resV)
                    m_pMesh->setNodeFlag(n, 0);

                if(resT)
                    m_pMesh->setNodeFlag(n, 1);

                if(resQ)
                    m_pMesh->setNodeFlag(n, 2);

                if(resT && resQ)
                    throw std::runtime_error("the boundary " + m_pMesh->getNodeType(n) +
                                             "has a BC for both T and Q. This is forbidden!");
            }
        }

        //Solve the Heat equation before or after the momentum continuity equation ?
        m_solveHeatFirst = m_solverParams[0].checkAndGet<bool>("solveHeatFirst");

        m_solveFunc = std::bind(&SolverIncompNewton::m_solveBoussinesq, this);
    }
    else if(m_pProblem->getID() == "Conduction")
    {
        //The heat equation
        m_pEquations.resize(1);

        bcFlags = {1, 2}; //Dirichlet and Neumann
        statesIndex = {0};
        m_pEquations[0] = std::make_unique<HeatEqIncompNewton>(
            m_pProblem, this, m_pMesh, m_solverParams, materialParams,
            bcFlags, statesIndex
        );

        //Set the right node flag if the boundary condition is present
        SolTable bcParamHeat = m_pEquations[0]->getBCParam(0);
        for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
        {
            const Node& node = m_pMesh->getNode(n);
            if(node.isBound())
            {
                bool resT = checkBC(bcParamHeat, n, node, "T", 1);
                bool resQ = checkBC(bcParamHeat, n, node, "Q", 1);

                if(resT)
                    m_pMesh->setNodeFlag(n, 1);

                if(resQ)
                    m_pMesh->setNodeFlag(n, 2);

                if(resT && resQ)
                    throw std::runtime_error("the boundary " + m_pMesh->getNodeType(n) +
                                             "has a BC for both T and Q. This is forbidden!");
            }
        }

        m_solveFunc = std::bind(&SolverIncompNewton::m_solveConduction, this);
    }

    //Loading time step parametrs
    m_adaptDT = m_solverParams[0].checkAndGet<bool>("adaptDT");
    m_maxDT = m_solverParams[0].checkAndGet<double>("maxDT");
    m_initialDT = m_solverParams[0].checkAndGet<double>("initialDT");
    m_coeffDTDecrease = m_solverParams[0].checkAndGet<double>("coeffDTDecrease");
    m_coeffDTincrease = m_solverParams[0].checkAndGet<double>("coeffDTincrease");

    m_timeStep = m_initialDT;

    //Should we compute the normals and the curvature ?
    m_pMesh->setComputeNormalCurvature(false);
    for(auto& pEq : m_pEquations)
    {
        if(pEq->isNormalCurvNeeded())
        {
            m_pMesh->setComputeNormalCurvature(true);
            break;
        }
    }
}

SolverIncompNewton::~SolverIncompNewton()
{

}

void SolverIncompNewton::displayParams() const
{
    std::cout << "Increase dt factor: " << m_coeffDTincrease << "\n"
              << "Decrease dt factor: " << m_coeffDTDecrease << "\n"
              << "Maximum dt: " << m_maxDT << "\n"
              << "Initial dt: " << m_initialDT << std::endl;

    for(auto& pEquation : m_pEquations)
        pEquation->displayParams();
}

bool SolverIncompNewton::solveOneTimeStep()
{
    return m_solveFunc();
}

void SolverIncompNewton::computeNextDT()
{
    if(!m_solveSucceed)
    {
        if(m_adaptDT)
            m_timeStep /= m_coeffDTDecrease;
        else
            throw std::runtime_error("solving the precedent time step was not successful!");
    }
    else
    {
        if(m_adaptDT)
            m_timeStep = std::min(m_maxDT, m_timeStep*m_coeffDTincrease);
    }
}

bool SolverIncompNewton::m_solveIncompNewtonNoT()
{
    m_solveSucceed = m_pEquations[0]->solve();

    if(m_solveSucceed)
    {
        m_pProblem->updateTime(m_timeStep);
        m_pMesh->remesh(m_pProblem->isOutputVerbose());
    }

    return m_solveSucceed;
}

bool SolverIncompNewton::m_solveBoussinesq()
{
    if(m_solveHeatFirst && !m_pEquations[1]->solve())
    {
        m_solveSucceed = false;
        return m_solveSucceed;
    }

    m_solveSucceed = m_pEquations[0]->solve();
    if(!m_solveSucceed)
        return m_solveSucceed;

    if(!m_solveHeatFirst && !m_pEquations[1]->solve())
    {
        m_solveSucceed = false;
        return m_solveSucceed;
    }

    if(m_solveSucceed)
    {
        m_pProblem->updateTime(m_timeStep);
        m_pMesh->remesh(m_pProblem->isOutputVerbose());
    }

    return m_solveSucceed;
}

bool SolverIncompNewton::m_solveConduction()
{
    m_solveSucceed = m_pEquations[0]->solve();

    if(m_solveSucceed)
        m_pProblem->updateTime(m_timeStep);

    return m_solveSucceed;
}
