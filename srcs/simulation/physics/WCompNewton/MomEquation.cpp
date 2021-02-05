#include "MomEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"
#include "../../utility/Clock.hpp"

MomEqWCompNewton::MomEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                                     std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                                     const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex) :
Equation(pProblem, pSolver, pMesh, solverParams, materialParams, bcFlags, statesIndex, "MomEq")
{
    m_mu = m_materialParams[0].checkAndGet<double>("mu");
    m_gamma = m_materialParams[0].checkAndGet<double>("gamma");

    if(bcFlags.size() != 1)
        throw std::runtime_error("the " + getID() + " only required one flag for one possible boundary condition!");

    if(m_pProblem->getID() == "BoussinesqWC")
    {
        m_alpha = m_materialParams[0].checkAndGet<double>("alpha");
        m_Tr = m_materialParams[0].checkAndGet<double>("Tr");

        if(statesIndex.size() != 5)
            throw std::runtime_error("the " + getID() + " equation requires 5 statesIndex: beginning of (u,v,w), beginning of (ax, ay, az), p and rho and T");
    }
    else
    {
        if(statesIndex.size() != 4)
            throw std::runtime_error("the " + getID() + "  equation requires 4 statesIndex: beginning of (u,v,w), beginning of (ax, ay, az), p and rho");
    }

    Eigen::MatrixXd ddev;
    if(m_pMesh->getDim() == 2)
    {
        ddev.resize(3, 3);
        ddev <<  4.0/3, -2.0/3, 0,
                -2.0/3,  4.0/3, 0,
                     0,      0, 1;
    }
    else
    {
        ddev.resize(6, 6);
        ddev <<  4.0/3, -2.0/3, -2.0/3, 0, 0, 0,
                -2.0/3,  4.0/3, -2.0/3, 0, 0, 0,
                -2.0/3, -2.0/3,  4.0/3, 0, 0, 0,
                     0,      0,      0, 1, 0, 0,
                     0,      0,      0, 0, 1, 0,
                     0,      0,      0, 0, 0, 1;
    }
    m_pMatBuilder->setddev(ddev);

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

    m_pMatBuilder->setMcomputeFactor([&](const Element& element, const Eigen::MatrixXd& N) -> double {
        return (N*getElementState(m_pMesh, element, m_statesIndex[3])).value();
    });

    m_pMatBuilder->setKcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
        return m_mu;
    });

    m_pMatBuilder->setDcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
        return 1;
    });

    if(m_pProblem->getID() == "BoussinesqWC")
    {
        m_pMatBuilder->setFcomputeFactor([&](const Element& element, const Eigen::MatrixXd&  N, const Eigen::MatrixXd& /** B **/) -> double {
            double rho = (N*getElementState(m_pMesh, element, m_statesIndex[3])).value();
            double T = (N*getElementState(m_pMesh, element, m_statesIndex[4])).value();
            return rho*(1 - m_alpha*(T - m_Tr));
        });
    }
    else
    {
        m_pMatBuilder->setFcomputeFactor([&](const Element& element, const Eigen::MatrixXd&  N, const Eigen::MatrixXd& /** B **/) -> double {
            return (N*getElementState(m_pMesh, element, m_statesIndex[3])).value();
        });
    }

    auto bodyForce = m_equationParams[0].checkAndGet<std::vector<double>>("bodyForce");
    if(bodyForce.size() != m_pMesh->getDim())
        throw std::runtime_error("the body force vector has not the right dimension!");

    m_bodyForce = Eigen::Map<Eigen::VectorXd>(bodyForce.data(), bodyForce.size());

    m_needNormalCurv = (m_gamma < 1e-15) ? false : true;
}

MomEqWCompNewton::~MomEqWCompNewton()
{

}

void MomEqWCompNewton::displayParams() const
{
    std::cout << "Momentum equation parameters:\n"
              << " * Viscosity: " << m_mu << " Pa s\n"
              << " * Surface Tension: " << m_gamma << " Nm" << std::endl;

    if(m_pMesh->getDim() == 2)
        std::cout << " * Body force: (" << m_bodyForce[0] << ", " << m_bodyForce[1] << ")" << std::endl;
    else
        std::cout << " * Body force: (" << m_bodyForce[0] << ", " << m_bodyForce[1] << "," << m_bodyForce[2] << ")" << std::endl;
}

double MomEqWCompNewton::getSpeedEquiv(double /** he **/, const Node& node)
{
    double nodeU = 0;
    for(unsigned int d = 0 ; d < m_pMesh->getDim() ; ++d)
    {
        nodeU += node.getState(d)*node.getState(d);
    }
    nodeU = std::sqrt(nodeU);

    return nodeU;
}

bool MomEqWCompNewton::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Momentum Equation" << std::endl;

    Eigen::VectorXd qV1half = getQFromNodesStates(m_pMesh, m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim() - 1);
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> invM; //The mass matrix for momentum equation.
    Eigen::VectorXd F;                                  //The rhs of the momentum equation.

    m_buildSystem(invM, F);
    m_applyBC(invM, F);
    Eigen::VectorXd qAcc = invM*F;

    Eigen::VectorXd qV = qV1half + 0.5*m_pSolver->getTimeStep()*qAcc;
    setNodesStatesfromQ(m_pMesh, qV, m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim() - 1);
    setNodesStatesfromQ(m_pMesh, qAcc, m_statesIndex[1], m_statesIndex[1] + m_pMesh->getDim() - 1);

    return true;
}

void MomEqWCompNewton::m_buildSystem(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F)
{
    const unsigned short dim = m_pMesh->getDim();
    const std::size_t elementsCount = m_pMesh->getElementsCount();
    const std::size_t nodesCount = m_pMesh->getNodesCount();
    const std::size_t noPerEl = m_pMesh->getNodesPerElm();

    invM.resize(dim*nodesCount); invM.setZero();
    std::vector<Eigen::MatrixXd> Me(elementsCount);
    for(auto& mat: Me)
        mat.resize(dim*noPerEl, dim*noPerEl);

    F.resize(dim*nodesCount); F.setZero();
    std::vector<Eigen::VectorXd> FTote(elementsCount);
    for(auto& vec: FTote)
        vec.resize(dim*noPerEl);

    Clock aClock;
    aClock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < elementsCount ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::VectorXd Rho = getElementState(m_pMesh, element, m_statesIndex[3]);

        Eigen::VectorXd V(noPerEl*dim);
        if(dim == 2)
        {
            V << getElementState(m_pMesh, element, m_statesIndex[0]),
                 getElementState(m_pMesh, element, m_statesIndex[0] + 1);
        }
        else
        {
             V << getElementState(m_pMesh, element, m_statesIndex[0]),
                  getElementState(m_pMesh, element, m_statesIndex[0] + 1),
                  getElementState(m_pMesh, element, m_statesIndex[0] + 2);
        }

        Eigen::VectorXd P = getElementState(m_pMesh, element, m_statesIndex[2]);

        Eigen::MatrixXd gradNe = m_pMatBuilder->getGradN(element);
        Eigen::MatrixXd Be = m_pMatBuilder->getB(gradNe);

        Eigen::MatrixXd MeTemp = m_pMatBuilder->getM(element);
        Me[elm] = MatrixBuilder::diagBlock(MeTemp, dim);
        MatrixBuilder::lump(Me[elm]);
        Eigen::MatrixXd Ke = m_pMatBuilder->getK(element, Be);
        Eigen::MatrixXd De = m_pMatBuilder->getD(element, Be);
        Eigen::VectorXd Fe = m_pMatBuilder->getF(element, m_bodyForce, Be);

        FTote[elm] = -Ke*V + De.transpose()*P + Fe;
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());

    auto& invMDiag = invM.diagonal();

    for(std::size_t elm = 0 ; elm < elementsCount ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                /********************************************************************
                                             Build M
                ********************************************************************/
                invMDiag[element.getNodeIndex(i) + d*nodesCount] += Me[elm](i + d*noPerEl, i + d*noPerEl);

                /************************************************************************
                                                Build f
                ************************************************************************/
                F(element.getNodeIndex(i) + d*nodesCount) += FTote[elm](i + d*noPerEl);
            }
        }
    }

    MatrixBuilder::inverse(invM);
}

void MomEqWCompNewton::m_applyBC(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F)
{
    assert(m_pMesh->getNodesCount() != 0);

    const uint8_t dim = m_pMesh->getDim();
    const std::size_t nodesCount = m_pMesh->getNodesCount();

    auto& invMDiag = invM.diagonal();

    //Do not parallelize this
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);

        if(node.isFree() && !node.isBound())
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                F(n + d*nodesCount) = m_bodyForce[d];
                invMDiag[n + d*nodesCount] = 1;
            }
        }
        else if(node.isBound())
        {
            if(node.getFlag(m_bcFlags[0]))
            {
                std::vector<double> result;
                result = m_bcParams[0].call<std::vector<double>>(m_pMesh->getNodeType(n) + "V",
                                                        node.getPosition(),
                                                        m_pMesh->getBoundNodeInitPos(n),
                                                        m_pProblem->getCurrentSimTime() +
                                                        m_pSolver->getTimeStep());

                for(unsigned short d = 0 ; d < dim ; ++d)
                {
                    F(n + d*nodesCount) = result[d];
                    invMDiag[n + d*nodesCount] = 1;
                }
            }
        }
    }
}
