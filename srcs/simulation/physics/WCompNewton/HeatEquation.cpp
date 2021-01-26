#include "HeatEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

HeatEqWCompNewton::HeatEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                                     std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                                     const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex) :
Equation(pProblem, pSolver, pMesh, solverParams, materialParams, bcFlags, statesIndex, "HeatEq")
{
    m_k = m_materialParams[0].checkAndGet<double>("k");
    m_cv = m_materialParams[0].checkAndGet<double>("cv");

    if(bcFlags.size() != 2)
        throw std::runtime_error("the " + getID() + " equation require two flags for two possible boundary conditions!");

    if(statesIndex.size() != 2)
        throw std::runtime_error("the " + getID() + " equation require two states index describing the T and rho state !");


    m_pMatBuilder->setMcomputeFactor([&](const Element& element, const Eigen::MatrixXd& N) -> double {
        double rho = (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
        return m_cv*rho;
    });

    m_pMatBuilder->setLcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
        return m_k;
    });

    m_needNormalCurv = false;
}

HeatEqWCompNewton::~HeatEqWCompNewton()
{

}

void HeatEqWCompNewton::displayParams() const
{
    std::cout << "Heat equation parameters:\n"
              << " * Specific heat capacity: " << m_cv << " J/(kg K)\n"
              << " * Heat conduction: " << m_k << " W/(mK)" << std::endl;
}

double HeatEqWCompNewton::getSpeedEquiv(double he, const Node& /** node **/)
{
    return 2*m_k/he;
}

bool HeatEqWCompNewton::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Heat Equation" << std::endl;

    Eigen::DiagonalMatrix<double,Eigen::Dynamic> invM; //The mass matrix of the continuity.
    Eigen::VectorXd F; //The mass matrix of the continuity.
    m_buildSystem(invM, F);
    m_applyBC(invM, F);

    Eigen::VectorXd qT = invM*F;

    setNodesStatesfromQ(m_pMesh, qT, m_statesIndex[0], m_statesIndex[0]);

    return true;
}

void HeatEqWCompNewton::m_applyBC(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F)
{
    auto& invMDiag = invM.diagonal();

    const std::size_t nodesCount = m_pMesh->getNodesCount();

    const std::size_t facetsCount = m_pMesh->getFacetsCount();
    const std::size_t noPerFacet = m_pMesh->getNodesPerFacet();

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree())
        {
            F(n) = node.getState(m_statesIndex[0]);
            invMDiag[n] = 1;
        }
        else if(node.getFlag(m_bcFlags[0]))
        {
            std::vector<double> result;
            result = m_bcParams[0].call<std::vector<double>>(m_pMesh->getNodeType(n) + "T",
                                                             node.getPosition(),
                                                             m_pMesh->getBoundNodeInitPos(n),
                                                             m_pProblem->getCurrentSimTime() +
                                                             m_pSolver->getTimeStep());
            F(n) = result[0];
            invMDiag[n] = 1;
        }
    }
}

void HeatEqWCompNewton::m_buildSystem(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F)
{
    const unsigned short dim = m_pMesh->getDim();
    const std::size_t elementsCount = m_pMesh->getElementsCount();
    const std::size_t nodesCount = m_pMesh->getNodesCount();

    invM.resize(nodesCount); invM.setZero();
    std::vector<Eigen::MatrixXd> Me(elementsCount);

    F.resize(nodesCount); F.setZero();
    std::vector<Eigen::VectorXd> FTote(elementsCount);

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < elementsCount ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::VectorXd T = getElementState(m_pMesh, element, m_statesIndex[0]);

        Eigen::MatrixXd gradNe = m_pMatBuilder->getGradN(element);
        Eigen::MatrixXd Be = m_pMatBuilder->getB(gradNe);
        Me[elm] = m_pMatBuilder->getM(element);
        MatrixBuilder::lump(Me[elm]);
        Eigen::MatrixXd Le = m_pMatBuilder->getL(element, Be, gradNe);

        FTote[elm] = - m_pSolver->getTimeStep()*Le*T + Me[elm]*T;
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());

    auto& invMDiag = invM.diagonal();

    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
        {
            invMDiag[element.getNodeIndex(i)] += Me[elm].diagonal()[i];

            F(element.getNodeIndex(i)) += FTote[elm](i);
        }
    }

    MatrixBuilder::inverse(invM);
}
