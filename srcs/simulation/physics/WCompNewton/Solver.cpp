#include "Problem.hpp"
#include "Solver.hpp"
#include "MomEquation.hpp"
#include "ContEquation.hpp"
#include "HeatEquation.hpp"
#include "../../utility/StatesFromToQ.hpp"
#include "../../utility/Clock.hpp"


SolverWCompNewton::SolverWCompNewton(Problem* pProblem, Mesh* pMesh, std::vector<SolTable> problemParams):
Solver(pProblem, pMesh, problemParams)
{
	//Check if the asked problem and solver are supported
    if(m_pProblem->getID() != "WCompNewtonNoT" && m_pProblem->getID() != "BoussinesqWC")
        throw std::runtime_error("this solver cannot be used with problem whose id is " + m_pProblem->getID());

    if(m_id != "CDS")
        throw std::runtime_error("this solver does not know id " + m_id);

    //Load material params for equations
    std::vector<SolTable> materialParams(problemParams.size());
    for(std::size_t i = 0 ; i < problemParams.size() ; ++i)
    {
        materialParams[i] = SolTable("Material", problemParams[i]);
    }

    unsigned int dim = m_pMesh->getDim();

	//Load equations depending of the problem
    std::vector<unsigned short> bcFlags;
    std::vector<unsigned int> statesIndex;
    if(m_pProblem->getID() == "WCompNewtonNoT")
    {
        //Momentum and continuity equations
        m_pEquations.resize(2);

		bcFlags = {0};
		statesIndex = {dim, dim + 1, 0};
		m_pEquations[0] = std::make_unique<ContEqWCompNewton>(
			m_pProblem, this, m_pMesh, m_solverParams, materialParams,
			bcFlags, statesIndex
		);

		bcFlags = {0};
		statesIndex = {0, dim + 2, dim, dim + 1};
		m_pEquations[1] = std::make_unique<MomEqWCompNewton>(
			m_pProblem, this, m_pMesh, m_solverParams, materialParams,
			bcFlags, statesIndex
		);

        //Set the right node flag if the boundary condition is present
        SolTable bcParam = m_pEquations[1]->getBCParam(0);
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

        m_solveFunc = std::bind(&SolverWCompNewton::m_solveWCompNewtonNoT, this);
    }
    else if(m_pProblem->getID() == "BoussinesqWC")
    {
        //The momentum, continuity and the heat equation
        m_pEquations.resize(3);

        bcFlags = {0};
		statesIndex = {dim, dim + 1, 0};
		m_pEquations[0] = std::make_unique<ContEqWCompNewton>(
			m_pProblem, this, m_pMesh, m_solverParams, materialParams,
			bcFlags, statesIndex
		);

		bcFlags = {0};
		statesIndex = {0, dim + 2, dim, dim + 1, 2*dim + 2};
		m_pEquations[1] = std::make_unique<MomEqWCompNewton>(
			m_pProblem, this, m_pMesh, m_solverParams, materialParams,
			bcFlags, statesIndex
		);

        bcFlags = {1, 2}; //Dirichlet and Neumann
        statesIndex = {2*dim + 2, dim + 1};
        m_pEquations[2] = std::make_unique<HeatEqWCompNewton>(
            m_pProblem, this, m_pMesh, m_solverParams, materialParams,
            bcFlags, statesIndex
        );

        //Set the right node flag if the boundary condition is present
        SolTable bcParamMomCont = m_pEquations[1]->getBCParam(0);
        SolTable bcParamHeat = m_pEquations[2]->getBCParam(0);
        for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
        {
            const Node& node = m_pMesh->getNode(n);
            if(node.isBound())
            {
                bool resV = checkBC(bcParamMomCont, n, node, "V", m_pMesh->getDim());
                bool resT = checkBC(bcParamHeat, n, node, "T", m_pMesh->getDim());
                bool resQ = checkBC(bcParamHeat, n, node, "Q", m_pMesh->getDim());

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

        m_solveFunc = std::bind(&SolverWCompNewton::m_solveBoussinesqWC, this);
	}

    m_adaptDT = m_solverParams[0].checkAndGet<bool>("adaptDT");
    m_maxDT = m_solverParams[0].checkAndGet<double>("maxDT");
    m_initialDT = m_solverParams[0].checkAndGet<double>("initialDT");
    m_securityCoeff = m_solverParams[0].checkAndGet<double>("securityCoeff");

    m_nextTimeToRemesh = m_maxDT;

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

SolverWCompNewton::~SolverWCompNewton()
{

}

void SolverWCompNewton::displayParams() const
{
    std::cout << "Maximum dt: " << m_maxDT << "\n"
              << "Initial dt: " << m_initialDT << "\n"
              << "Security coeff: " << m_securityCoeff << std::endl;

    for(auto& pEquation : m_pEquations)
        pEquation->displayParams();
}

bool SolverWCompNewton::solveOneTimeStep()
{
    return m_solveFunc();
}

void SolverWCompNewton::computeNextDT()
{
    if(m_adaptDT)
    {
        m_timeStep = std::numeric_limits<double>::max();
        for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
        {
            const Element& element = m_pMesh->getElement(elm);
            double he = 2*element.getRin();

            double maxSpeed = 0;
            for(std::size_t n = 0 ; n < m_pMesh->getNodesPerElm() ; ++n)
            {
                const Node& node = element.getNode(n);

                double c = m_pEquations[0]->getSpeedEquiv(he, node);
                double u = m_pEquations[1]->getSpeedEquiv(he, node);

                maxSpeed = std::max(std::max(u, c), maxSpeed);
            }
            m_timeStep = std::min(m_timeStep, m_securityCoeff*he/maxSpeed);
        }

        m_timeStep = std::min(m_timeStep, m_maxDT);
        if(std::isnan(m_timeStep) || std::isinf(m_timeStep))
            m_timeStep = m_maxDT;
    }
}

bool SolverWCompNewton::m_solveWCompNewtonNoT()
{
    unsigned int dim = m_pMesh->getDim();
    Eigen::VectorXd qVPrev = getQFromNodesStates(m_pMesh, 0, dim - 1);           //The precedent speed.
    Eigen::VectorXd qAccPrev = getQFromNodesStates(m_pMesh, dim + 2, 2*dim + 1);  //The precedent acceleration.

    m_pEquations[0]->preCompute();

    Eigen::VectorXd qV1half = qVPrev + 0.5*m_timeStep*qAccPrev;

    setNodesStatesfromQ(m_pMesh, qV1half, 0, dim - 1);
    Eigen::VectorXd deltaPos = qV1half*m_timeStep;
    m_pMesh->updateNodesPosition(std::vector<double> (deltaPos.data(), deltaPos.data() + deltaPos.cols()*deltaPos.rows()));

    m_pEquations[0]->solve();
    m_pEquations[1]->solve();

    m_pProblem->updateTime(m_timeStep);
    if(m_pProblem->getCurrentSimTime() > m_nextTimeToRemesh)
    {
        m_pMesh->remesh(m_pProblem->isOutputVerbose());
        m_nextTimeToRemesh += m_maxDT;
    }

    return true;
}

bool SolverWCompNewton::m_solveBoussinesqWC()
{
    unsigned int dim = m_pMesh->getDim();
    Eigen::VectorXd qVPrev = getQFromNodesStates(m_pMesh, 0, dim - 1);           //The precedent speed.
    Eigen::VectorXd qAccPrev = getQFromNodesStates(m_pMesh, dim + 2, 2*dim + 1);  //The precedent acceleration.

    m_pEquations[0]->preCompute();

    Eigen::VectorXd qV1half = qVPrev + 0.5*m_timeStep*qAccPrev;

    setNodesStatesfromQ(m_pMesh, qV1half, 0, dim - 1);
    Eigen::VectorXd deltaPos = qV1half*m_timeStep;
    m_pMesh->updateNodesPosition(std::vector<double> (deltaPos.data(), deltaPos.data() + deltaPos.cols()*deltaPos.rows()));

    m_pEquations[2]->solve();
    m_pEquations[0]->solve();
    m_pEquations[1]->solve();

    m_pProblem->updateTime(m_timeStep);
    if(m_pProblem->getCurrentSimTime() > m_nextTimeToRemesh)
    {
        m_pMesh->remesh(m_pProblem->isOutputVerbose());
        m_nextTimeToRemesh += m_maxDT;
    }

    return true;
}
