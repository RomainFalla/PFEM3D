#include "HeatEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

HeatEqIncompNewton::HeatEqIncompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                                     std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                                     const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex) :
Equation(pProblem, pSolver, pMesh, solverParams, materialParams, bcFlags, statesIndex, "HeatEq")
{
    m_rho = m_materialParams[0].checkAndGet<double>("rho");
    m_k = m_materialParams[0].checkAndGet<double>("k");
    m_cv = m_materialParams[0].checkAndGet<double>("cv");

    if(bcFlags.size() != 2)
        throw std::runtime_error("the " + getID() + " equation require two flags for two possible boundary conditiosn!");

    if(statesIndex.size() != 1)
        throw std::runtime_error("the " + getID() + " equation require one state index describing the T state !");


    m_pMatBuilder->setMcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/) -> double {
        return m_rho*m_cv;
    });

    m_pMatBuilder->setLcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
        return m_k;
    });

    m_pMatBuilder->setQFunc([&](const Facet& facet, const std::array<double, 3>& gp) -> Eigen::VectorXd {
        std::array<double, 3> pos = facet.getPosFromGP(gp);
        std::vector<double> result;
            result = m_bcParams[0].call<std::vector<double>>(m_pMesh->getNodeType(facet.getNodeIndex(0)) + "Q",
                                                             pos,
                                                             pos, //TO DO: fix
                                                             m_pProblem->getCurrentSimTime() +
                                                             m_pSolver->getTimeStep());

        return Eigen::VectorXd::Map(result.data(), result.size());
    });

    unsigned int maxIter = m_equationParams[0].checkAndGet<unsigned int>("maxIter");
    double minRes = m_equationParams[0].checkAndGet<double>("minRes");

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

        for(auto i = 0 ; i < qIter.rows() ; ++i)
        {
            num += (qIter(i) - qIterPrev(i))*(qIter(i) - qIterPrev(i));
            den += qIterPrev(i)*qIterPrev(i);
        }

        return std::sqrt(num/den);
    }, maxIter, minRes);

    m_needNormalCurv = true;
}

HeatEqIncompNewton::~HeatEqIncompNewton()
{

}

void HeatEqIncompNewton::displayParams() const
{
    std::cout << "Heat equation parameters:\n"
              << " * Specific heat capacity: " << m_cv << " J/(kg K)\n"
              << " * Heat conduction: " << m_k << " W/(mK)" << std::endl;

    m_pPicardAlgo->displayParams();
}

bool HeatEqIncompNewton::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Heat Equation" << std::endl;

    Eigen::VectorXd qPrev = getQFromNodesStates(m_pMesh, m_statesIndex[0], m_statesIndex[0]);
    return m_pPicardAlgo->solve(m_pMesh, qPrev, m_pProblem->isOutputVerbose());
}

void HeatEqIncompNewton::m_buildAb(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::VectorXd& qPrev)
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

        Eigen::MatrixXd gradNe = m_pMatBuilder->getGradN(element);
        Eigen::MatrixXd Be = m_pMatBuilder->getB(gradNe);
        Eigen::MatrixXd Me_dt = (1/dt)*m_pMatBuilder->getM(element);
        Eigen::MatrixXd Le = m_pMatBuilder->getL(element, Be, gradNe);

        Eigen::VectorXd thetaPrev(noPerEl);
         if(dim == 2)
            thetaPrev << qPrev[element.getNodeIndex(0)], qPrev[element.getNodeIndex(1)], qPrev[element.getNodeIndex(2)];
        else
            thetaPrev << qPrev[element.getNodeIndex(0)], qPrev[element.getNodeIndex(1)], qPrev[element.getNodeIndex(2)], qPrev[element.getNodeIndex(3)];

        Eigen::VectorXd MThetaPreve_dt = Me_dt*thetaPrev;

        std::size_t countA = 0;
        std::size_t countb = 0;

         for(unsigned short i = 0 ; i < noPerEl ; ++i)
        {
            const Node& ni = element.getNode(i);

            for(unsigned short j = 0 ; j < noPerEl ; ++j)
            {
                if(!ni.getFlag(m_bcFlags[0]))
                {
                    if(!ni.isFree())
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(element.getNodeIndex(i),
                                                   element.getNodeIndex(j),
                                                   Me_dt(i, j));
                    }
                    countA++;

                    if(!ni.isFree())
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(element.getNodeIndex(i),
                                                   element.getNodeIndex(j),
                                                   Le(i, j));
                    }
                    countA++;
                }
            }

            /************************************************************************
                                            Build h
            ************************************************************************/
            indexb[noPerEl*elm + countb] = std::make_pair(element.getNodeIndex(i), MThetaPreve_dt[i]);

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

        if(node.getFlag(m_bcFlags[0]) || node.isFree())
        {
            indexA.push_back(Eigen::Triplet<double>(n, n, 1));
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

void HeatEqIncompNewton::m_applyBC(Eigen::VectorXd& b, const Eigen::VectorXd& qPrev)
{
    const std::size_t nodesCount = m_pMesh->getNodesCount();

    const std::size_t facetsCount = m_pMesh->getFacetsCount();
    const std::size_t noPerFacet = m_pMesh->getNodesPerFacet();

    for(std::size_t f = 0 ; f < facetsCount ; ++f)
    {
        const Facet& facet = m_pMesh->getFacet(f);
        bool boundaryQ = true;
        for(uint8_t n = 0 ; n < noPerFacet ; ++n)
        {
            if(!m_pMesh->getNode(facet.getNodeIndex(n)).getFlag(m_bcFlags[1]))
            {
                boundaryQ = false;
                break;
            }
        }
        if(!boundaryQ)
            continue;

        Eigen::VectorXd qn = m_pMatBuilder->getQN(facet);

        for(unsigned short i = 0 ; i < noPerFacet ; ++i)
        {
            b(facet.getNodeIndex(i)) -= qn[i]; //See Galerkin formulation
        }
    }

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree())
        {
            b(n) = qPrev[n];
        }
        else if(node.getFlag(m_bcFlags[0]))
        {
            std::vector<double> result;
            result = m_bcParams[0].call<std::vector<double>>(m_pMesh->getNodeType(n) + "T",
                                                             node.getPosition(),
                                                             m_pMesh->getBoundNodeInitPos(n),
                                                             m_pProblem->getCurrentSimTime() +
                                                             m_pSolver->getTimeStep());
            b(n) = result[0];
        }
    }
}

void HeatEqIncompNewton::m_executeTask(const Eigen::VectorXd& qIter)
{
    setNodesStatesfromQ(m_pMesh, qIter, m_statesIndex[0], m_statesIndex[0]);
}
