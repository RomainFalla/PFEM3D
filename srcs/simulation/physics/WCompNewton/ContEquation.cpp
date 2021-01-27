#include "ContEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

ContEqWCompNewton::ContEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                                     std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                                     const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex) :
Equation(pProblem, pSolver, pMesh, solverParams, materialParams, bcFlags, statesIndex, "ContEq")
{
    m_K0 = m_materialParams[0].checkAndGet<double>("K0");
    m_K0p = m_materialParams[0].checkAndGet<double>("K0p");
    m_rhoStar = m_materialParams[0].checkAndGet<double>("rhoStar");
    m_strongContinuity = m_equationParams[0].checkAndGet<bool>("strongContinuity");

    Eigen::VectorXd m;
    if(m_pMesh->getDim() == 2)
    {
        m.resize(3);
        m << 1, 1, 0;
    }
    else
    {
        m.resize(6);
        m << 1, 1, 1, 0, 0, 0;
    }
    m_pMatBuilder->setm(m);

    if(m_statesIndex .size()!= 3)
        throw std::runtime_error("the " + getID() + " equation requires 3 statesIndex: one index for the p unknown, one for the rho unknown, and one for the beginning of the (u,v,w) unknown!");

    m_pMatBuilder->setMcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/) -> double {
        return 1;
    });

    m_pMatBuilder->setDcomputeFactor([&](const Element& element, const Eigen::MatrixXd& N, const Eigen::MatrixXd& /** B **/) -> double {
        return (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
    });

    m_needNormalCurv = false;
}

ContEqWCompNewton::~ContEqWCompNewton()
{

}

void ContEqWCompNewton::displayParams() const
{
    std::cout << "Continuity equation parameters:\n"
              << " * K0: " << m_K0 << " Pa\n"
              << " * K0': " << m_K0p << " \n"
              << " * Density at zero p: " << m_rhoStar << " Kg/m^3\n"
              << " * Strong continuity: " << std::boolalpha << m_strongContinuity << std::endl;
}

double ContEqWCompNewton::getSpeedEquiv(double /** he **/, const Node& node)
{
    return (m_K0 + m_K0p*node.getState(m_pMesh->getDim()))/node.getState(m_pMesh->getDim() + 1);
}

bool ContEqWCompNewton::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Continuity Equation" << std::endl;

    Eigen::DiagonalMatrix<double,Eigen::Dynamic> invM; //The mass matrix of the continuity.
    m_buildSystem(invM, m_F0);
    m_applyBC(invM, m_F0);

    Eigen::VectorXd qRho(m_pMesh->getNodesCount());
    qRho = invM*m_F0;
    setNodesStatesfromQ(m_pMesh, qRho, m_statesIndex[1], m_statesIndex[1]);

    Eigen::VectorXd qP = m_getPFromRhoTaitMurnagham(qRho);
    setNodesStatesfromQ(m_pMesh, qP, m_statesIndex[0], m_statesIndex[0]);

    return true;
}

void ContEqWCompNewton::preCompute()
{
    if(m_strongContinuity)
        m_buildF0(m_F0);
}

void ContEqWCompNewton::m_applyBC(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F0)
{
    auto& invMDiag = invM.diagonal();

    for (std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);

        if(node.isFree())
        {
            F0(n) = m_rhoStar;

            invMDiag[n] = 1;
        }
    }
}

void ContEqWCompNewton::m_buildF0(Eigen::VectorXd& F0)
{
    F0.resize(m_pMesh->getNodesCount()); F0.setZero();

    std::vector<Eigen::VectorXd> F0e(m_pMesh->getElementsCount());

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::VectorXd Rho = getElementState(m_pMesh, element, m_statesIndex[1]);

        Eigen::MatrixXd Mrhoe = m_pMatBuilder->getM(element);

        F0e[elm] = Mrhoe*Rho;
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());

    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(uint8_t i = 0 ; i < m_pMesh->getNodesPerElm() ; ++i)
            F0(element.getNodeIndex(i)) += F0e[elm](i);
    }
}

void ContEqWCompNewton::m_buildSystem(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F0)
{
    const unsigned short dim = m_pMesh->getDim();
    const unsigned short noPerEl = m_pMesh->getNodesPerElm();

    invM.resize(m_pMesh->getNodesCount()); invM.setZero();

    std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> MeLumped(m_pMesh->getElementsCount());
    std::vector<Eigen::VectorXd> F0e;

    if(!m_strongContinuity)
    {
        F0.resize(m_pMesh->getNodesCount()); F0.setZero();
        F0e.resize(m_pMesh->getElementsCount());
    }

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::MatrixXd Me = m_pMatBuilder->getM(element);
        MeLumped[elm] = MatrixBuilder::lump2(Me);

        if(!m_strongContinuity)
        {
            Eigen::VectorXd Rho = getElementState(m_pMesh, element, m_statesIndex[1]);
            Eigen::VectorXd V(noPerEl*dim);
            if(dim == 2)
            {
                V << getElementState(m_pMesh, element, m_statesIndex[2]),
                     getElementState(m_pMesh, element, m_statesIndex[2] + 1);
            }
            else
            {
                 V << getElementState(m_pMesh, element, m_statesIndex[2]),
                      getElementState(m_pMesh, element, m_statesIndex[2] + 1),
                      getElementState(m_pMesh, element, m_statesIndex[2] + 2);
            }

            Eigen::MatrixXd gradNe = m_pMatBuilder->getGradN(element);
            Eigen::MatrixXd Be = m_pMatBuilder->getB(gradNe);
            Eigen::MatrixXd Drhoe = m_pMatBuilder->getD(element, Be);

            F0e[elm] = MeLumped[elm]*Rho - m_pSolver->getTimeStep()*Drhoe*V;
        }
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());

    auto& invMDiag = invM.diagonal();

    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(unsigned short i = 0 ; i < noPerEl ; ++i)
        {
            invMDiag[element.getNodeIndex(i)] += MeLumped[elm].diagonal()[i];

            if(!m_strongContinuity)
            {
                F0(element.getNodeIndex(i)) += F0e[elm](i);
            }
        }
    }

    MatrixBuilder::inverse(invM);
}

Eigen::VectorXd ContEqWCompNewton::m_getPFromRhoTaitMurnagham(const Eigen::VectorXd& qRho)
{
    Eigen::VectorXd qP(m_pMesh->getNodesCount());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
    {
        double rho = qRho[n];
        qP[n] = (m_K0/m_K0p)*(std::pow(rho/m_rhoStar, m_K0p) - 1);
    }

    return qP;
}
