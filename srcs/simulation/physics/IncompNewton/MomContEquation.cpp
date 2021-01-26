#include "MomContEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

MomContEqIncompNewton::MomContEqIncompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                                     std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                                     const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex) :
Equation(pProblem, pSolver, pMesh, solverParams, materialParams, bcFlags, statesIndex, "MomContEq")
{
    m_rho = m_materialParams[0].checkAndGet<double>("rho");
    m_mu = m_materialParams[0].checkAndGet<double>("mu");
    m_gamma = m_materialParams[0].checkAndGet<double>("gamma");

    if(m_pProblem->getID() == "Boussinesq")
    {
        m_alpha = m_materialParams[0].checkAndGet<double>("alpha");
        m_Tr = m_materialParams[0].checkAndGet<double>("Tr");

        if(statesIndex.size() != 2)
            throw std::runtime_error("the " + getID() + " equation require two statesIndex describing the beginning of the states span and the temperature state!");
    }
    else
    {
        if(statesIndex.size() != 1)
            throw std::runtime_error("the " + getID() + " equation require one statesIndex describing the beginning of the states span!");
    }

    if(bcFlags.size() != 1)
        throw std::runtime_error("the " + getID() + " equation require one flag for one possible boundary condition!");



    Eigen::MatrixXd ddev;
    if(m_pMesh->getDim() == 2)
    {
        ddev.resize(3, 3);
        ddev << 2, 0, 0,
                0, 2, 0,
                0, 0, 1;
    }
    else
    {
        ddev.resize(6, 6);
        ddev << 2, 0, 0, 0, 0, 0,
                0, 2, 0, 0, 0, 0,
                0, 0, 2, 0, 0, 0,
                0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 1;
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

    m_pMatBuilder->setMcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/) -> double {
        return m_rho;
    });

    m_pMatBuilder->setMGammacomputeFactor([&](const Facet& /** facet **/, const Eigen::MatrixXd& /** N **/) -> double {
        return m_gamma;
    });

    m_pMatBuilder->setKcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
        return m_mu;
    });

    m_pMatBuilder->setDcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
        return 1;
    });

    m_pMatBuilder->setCcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
        return 1;
    });

    m_pMatBuilder->setLcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
        return 1/m_rho;
    });

    if(m_pProblem->getID() == "Boussinesq")
    {
        m_pMatBuilder->setFcomputeFactor([&](const Element& element, const Eigen::MatrixXd&  N, const Eigen::MatrixXd& /** B **/) -> double {
            double T = (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
            return m_rho*(1 - m_alpha*(T - m_Tr));
        });

        m_pMatBuilder->setHcomputeFactor([&](const Element& element, const Eigen::MatrixXd& N, const Eigen::MatrixXd& /** B **/) -> double {
            double T = (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
            return (1 - m_alpha*(T - m_Tr));
        });
    }
    else
    {
        m_pMatBuilder->setFcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
            return m_rho;
        });

        m_pMatBuilder->setHcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
            return 1;
        });
    }

    unsigned int maxIter = m_equationParams[0].checkAndGet<unsigned int>("maxIter");
    double minRes = m_equationParams[0].checkAndGet<double>("minRes");

    auto bodyForce = m_equationParams[0].checkAndGet<std::vector<double>>("bodyForce");
    if(bodyForce.size() != m_pMesh->getDim())
        throw std::runtime_error("the body force vector has not the right dimension!");

    m_bodyForce = Eigen::Map<Eigen::VectorXd>(bodyForce.data(), bodyForce.size());

    m_pPicardAlgo = std::make_unique<PicardAlgo>([this](auto& A, auto& b, const auto& qPrev){
        this->m_buildAb(A, b, qPrev);
    },
    [this](auto& b, const auto& qPrev){
        this->m_applyBC(b, qPrev);
    },
    [this](const auto& q){
        this->m_executeTask(q);
    },
    [this](const auto& qIter, const auto& qIterPrev) -> double {
        double num = 0, den = 0;
        Mesh* p_Mesh = this->m_pMesh;

        for(std::size_t n = 0 ; n < p_Mesh->getNodesCount() ; ++n)
        {
            const Node& node = p_Mesh->getNode(n);

            if(!node.isFree())
            {
                for(unsigned short d = 0 ; d < p_Mesh->getDim() ; ++d)
                {
                    num += (qIter(n + d*p_Mesh->getNodesCount()) - qIterPrev(n + d*p_Mesh->getNodesCount()))*(qIter(n + d*p_Mesh->getNodesCount()) - qIterPrev(n + d*p_Mesh->getNodesCount()));
                    den += qIterPrev(n + d*p_Mesh->getNodesCount())*qIterPrev(n + d*p_Mesh->getNodesCount());
                }
            }
        }

        return std::sqrt(num/den);
    }, maxIter, minRes);

    m_needNormalCurv = (m_gamma < 1e-15) ? false : true;
}

MomContEqIncompNewton::~MomContEqIncompNewton()
{

}

void MomContEqIncompNewton::displayParams() const
{
    std::cout << "Momentum Continuity equation parameters: \n"
              << " * Density: " << m_rho << " kg/m^3\n"
              << " * Viscosity: " << m_mu << " Pa s\n"
              << " * Surface Tension: " << m_gamma << " Nm" << std::endl;

    if(m_pProblem->getID() == "Boussinesq")
    {
        std::cout << " * Thermal expansion coefficient: " << m_alpha << " 1/K\n"
                  << " * Reference temperature: " << m_Tr << " K" << std::endl;
    }

    if(m_pMesh->getDim() == 2)
        std::cout << " * Body force: (" << m_bodyForce[0] << ", " << m_bodyForce[1] << ")" << std::endl;
    else
        std::cout << " * Body force: (" << m_bodyForce[0] << ", " << m_bodyForce[1] << "," << m_bodyForce[2] << ")" << std::endl;

    m_pPicardAlgo->displayParams();
}

bool MomContEqIncompNewton::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Momentum Continuity Equation" << std::endl;

    Eigen::VectorXd qPrev = getQFromNodesStates(m_pMesh, m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim());
    return m_pPicardAlgo->solve(m_pMesh, qPrev, m_pProblem->isOutputVerbose());
}

void MomContEqIncompNewton::m_buildAb(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::VectorXd& qPrev)
{
    unsigned int dim = m_pMesh->getDim();
    unsigned int noPerEl = m_pMesh->getNodesPerElm();
    unsigned int tripletPerElm = (dim*noPerEl*noPerEl + dim*noPerEl*dim*noPerEl + 3*noPerEl*dim*noPerEl + noPerEl*noPerEl);
    unsigned int doubletPerElm = (2*dim*noPerEl + 2*noPerEl);
    std::size_t nElm = m_pMesh->getElementsCount();
    std::size_t nNodes = m_pMesh->getNodesCount();
    double dt = m_pSolver->getTimeStep();

    std::vector<Eigen::Triplet<double>> indexA(tripletPerElm*nElm);
    std::vector<std::pair<std::size_t, double>> indexb(doubletPerElm*nElm); b.setZero();

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        double tau = m_computeTauPSPG(element);
        Eigen::MatrixXd gradNe = m_pMatBuilder->getGradN(element);
        Eigen::MatrixXd Be = m_pMatBuilder->getB(gradNe);
        Eigen::MatrixXd Me_dt = (1/dt)*m_pMatBuilder->getM(element);
        Me_dt = MatrixBuilder::diagBlock(Me_dt, dim);
        Eigen::MatrixXd Ke = m_pMatBuilder->getK(element, Be);
        Eigen::MatrixXd De = m_pMatBuilder->getD(element, Be);
        Eigen::MatrixXd Ce_dt = (tau/dt)*m_pMatBuilder->getC(element, Be, gradNe);
        Eigen::MatrixXd Le = tau*m_pMatBuilder->getL(element, Be, gradNe);
        Eigen::VectorXd Fe = m_pMatBuilder->getF(element, m_bodyForce, Be);
        Eigen::VectorXd He = tau*m_pMatBuilder->getH(element, m_bodyForce, Be, gradNe);

        Eigen::VectorXd vPrev(noPerEl*dim);
        if(dim == 2)
            vPrev << qPrev[element.getNodeIndex(0)], qPrev[element.getNodeIndex(1)], qPrev[element.getNodeIndex(2)],
                     qPrev[element.getNodeIndex(0) + nNodes], qPrev[element.getNodeIndex(1) + nNodes], qPrev[element.getNodeIndex(2) + nNodes];
        else
            vPrev << qPrev[element.getNodeIndex(0)], qPrev[element.getNodeIndex(1)], qPrev[element.getNodeIndex(2)], qPrev[element.getNodeIndex(3)],
                     qPrev[element.getNodeIndex(0) + nNodes], qPrev[element.getNodeIndex(1) + nNodes], qPrev[element.getNodeIndex(2) + nNodes], qPrev[element.getNodeIndex(3) + nNodes],
                     qPrev[element.getNodeIndex(0) + 2*nNodes], qPrev[element.getNodeIndex(1) + 2*nNodes], qPrev[element.getNodeIndex(2) + 2*nNodes], qPrev[element.getNodeIndex(3) + 2*nNodes];

        Eigen::VectorXd MvPreve_dt = Me_dt*vPrev;
        Eigen::VectorXd CvPreve_dt = Ce_dt*vPrev;

        std::size_t countA = 0;
        std::size_t countb = 0;

        for(unsigned short i = 0 ; i < noPerEl ; ++i)
        {
            const Node& ni = m_pMesh->getNode(element.getNodeIndex(i));

            for(unsigned short j = 0 ; j < noPerEl ; ++j)
            {
                const Node& nj = m_pMesh->getNode(element.getNodeIndex(j));

                for(unsigned short d = 0 ; d < dim ; ++d)
                {
                    /********************************************************************
                                                 Build M/dt
                    ********************************************************************/
                    if(!(ni.isBound() || ni.isFree()))
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(element.getNodeIndex(i) + d*nNodes,
                                                   element.getNodeIndex(j) + d*nNodes,
                                                   Me_dt(i + d*noPerEl, j + d*noPerEl));
                    }
                    countA++;

                    /********************************************************************
                                                  Build K
                    ********************************************************************/
                    for(unsigned short d2 = 0 ; d2 < dim ; ++d2)
                    {

                        if(!(ni.isBound() || ni.isFree()))
                        {
                            indexA[tripletPerElm*elm + countA] =
                                Eigen::Triplet<double>(element.getNodeIndex(i) + d*nNodes,
                                                       element.getNodeIndex(j) + d2*nNodes,
                                                       Ke(i + d*noPerEl, j + d2*noPerEl));
                        }
                        countA++;
                    }

                    /********************************************************************
                                                  Build D
                    ********************************************************************/
                    if(!ni.isFree())
                    {
                        //Part D of A
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(element.getNodeIndex(i) + dim*nNodes,
                                                   element.getNodeIndex(j) + d*nNodes,
                                                   De(i, j + d*noPerEl));
                    }
                    countA++;

                    //Part -D^T of A
                    if(!(nj.isBound() || nj.isFree()))
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(element.getNodeIndex(j) + d*nNodes,
                                                   element.getNodeIndex(i) + dim*nNodes,
                                                   -De(i, j + d*noPerEl));
                    }
                    countA++;

                    /********************************************************************
                                                Build C/dt
                    ********************************************************************/
                    if(!ni.isFree())
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(element.getNodeIndex(i) + dim*nNodes,
                                                   element.getNodeIndex(j) + d*nNodes,
                                                   Ce_dt(i, j + d*noPerEl));
                    }
                    countA++;
                }

                /********************************************************************
                                            Build L
                ********************************************************************/
                if(!ni.isFree())
                {
                    indexA[tripletPerElm*elm + countA] =
                        Eigen::Triplet<double>(element.getNodeIndex(i) + dim*nNodes,
                                               element.getNodeIndex(j) + dim*nNodes,
                                               Le(i,j));
                }
                countA++;
            }

            /************************************************************************
                                              Build f
            ************************************************************************/
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                indexb[doubletPerElm*elm + countb] = std::make_pair(element.getNodeIndex(i) + d*nNodes, Fe(i + d*noPerEl));
                countb++;

                indexb[doubletPerElm*elm + countb] = std::make_pair(element.getNodeIndex(i) + d*nNodes, MvPreve_dt(i + d*noPerEl));
                countb++;
            }

            /************************************************************************
                                            Build h
            ************************************************************************/

            indexb[doubletPerElm*elm + countb] = std::make_pair(element.getNodeIndex(i) + dim*nNodes, He(i));
            countb++;

            indexb[doubletPerElm*elm + countb] = std::make_pair(element.getNodeIndex(i) + dim*nNodes, CvPreve_dt(i));
            countb++;
        }
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());

    //Best would be to know the number of nodes in which case :/
    //This can still be fasten using OpenMP but will never be as good as using []
    //with preallocated memory
    for(std::size_t n = 0 ; n < nNodes ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);

        if(node.isFree())
        {
            indexA.push_back(Eigen::Triplet<double>(n + dim*nNodes,
                                                    n + dim*nNodes,
                                                    1));
        }

        if(node.isBound() || node.isFree())
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                indexA.push_back(Eigen::Triplet<double>(n + d*nNodes,
                                                        n + d*nNodes,
                                                        1));
            }
        }
    }

    /********************************************************************************
                                        Compute A and b
    ********************************************************************************/
    A.setFromTriplets(indexA.begin(), indexA.end());

    for(const auto& doublet : indexb)
    {
        //std::cout << doublet.first << ", " << doublet.second << std::endl;
        b[doublet.first] += doublet.second;
    }
}

void MomContEqIncompNewton::m_applyBC(Eigen::VectorXd& b, const Eigen::VectorXd& qPrev)
{
    const uint8_t dim = m_pMesh->getDim();
    const std::size_t nodesCount = m_pMesh->getNodesCount();

    const std::size_t facetsCount = m_pMesh->getFacetsCount();
    const std::size_t noPerFacet = m_pMesh->getNodesPerFacet();

    for(std::size_t f = 0 ; f < facetsCount ; ++f)
    {
        if(m_gamma < 1e-15)
            continue;

        const Facet& facet = m_pMesh->getFacet(f);

        bool onFS = true;
        for(unsigned short n = 0 ; n < noPerFacet ; ++n)
        {
            if(!facet.getNode(n).isOnFreeSurface())
            {
                onFS = false;
                break;
            }
        }
        if(!onFS)
            continue;

        Eigen::MatrixXd MGamma = m_pMatBuilder->getMGamma(facet);
        MGamma = MatrixBuilder::diagBlock(MGamma, m_pMesh->getDim());

        Eigen::VectorXd kappa_n(dim*noPerFacet);
        for(uint8_t n = 0 ; n < noPerFacet; ++n)
        {
            double curvature = m_pMesh->getFreeSurfaceCurvature(facet.getNodeIndex(n));
            std::array<double, 3> normal = m_pMesh->getBoundFSNormal(facet.getNodeIndex(n));

            for(uint8_t d = 0 ; d < dim ; ++d)
                kappa_n(n + d*noPerFacet) = curvature*normal[d];
        }

        Eigen::VectorXd Ff = MGamma*kappa_n;

        for(unsigned short i = 0 ; i < noPerFacet ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
                b(facet.getNodeIndex(i) + d*nodesCount) += Ff[d*noPerFacet + i];
        }
    }

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree())
        {
            b(n + dim*nodesCount) = 0;

            if(!node.isBound())
            {
                for(uint8_t d = 0 ; d < dim  ; ++d)
                {
                    b(n + d*nodesCount) = qPrev(n + d*nodesCount) + m_pSolver->getTimeStep()*m_bodyForce[d];
                }
            }
        }

        if(node.isBound())
        {
            if(node.getFlag(m_bcFlags[0]))
            {
                std::vector<double> result;
                result = m_bcParams[0].call<std::vector<double>>(m_pMesh->getNodeType(n) + "V",
                                                        node.getPosition(),
                                                        m_pMesh->getBoundNodeInitPos(n),
                                                        m_pProblem->getCurrentSimTime() +
                                                        m_pSolver->getTimeStep());

                for(uint8_t d = 0 ; d < dim ; ++d)
                {
                    b(n + d*nodesCount) = result[d];
                }
            }
        }
    }
}

void MomContEqIncompNewton::m_executeTask(const Eigen::VectorXd& qIter)
{
    setNodesStatesfromQ(m_pMesh, qIter, m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim());
    Eigen::VectorXd deltaPos = qIter*m_pSolver->getTimeStep();
    m_pMesh->updateNodesPositionFromSave(std::vector<double> (deltaPos.data(), deltaPos.data() + m_pMesh->getDim()*m_pMesh->getNodesCount()));
}

double MomContEqIncompNewton::m_computeTauPSPG(const Element& element) const
{
    const double h = std::sqrt(m_pMesh->getRefElementSize(m_pMesh->getDim())*element.getDetJ()/M_PI);

    double U = 0;
    for (unsigned short n = 0 ; n < m_pMesh->getNodesPerElm() ; ++n)
    {
        const Node& node = m_pMesh->getNode(element.getNodeIndex(n));

        double nodeU = 0;
        for (unsigned short d = 0 ; d < m_pMesh->getDim() ; ++d)
        {
            nodeU += node.getState(d)*node.getState(d);
        }
        U += std::sqrt(nodeU);
    }
    U /= (m_pMesh->getDim() + 1);

    return 1/std::sqrt((2/m_pSolver->getTimeStep())*(2/m_pSolver->getTimeStep()) + (2*U/h)*(2*U/h)
                        + 9*(4*m_mu/(h*h*m_rho))*(4*m_mu/(h*h*m_rho)));
}
