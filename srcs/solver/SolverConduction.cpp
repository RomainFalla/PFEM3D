#include <chrono>
#include <iomanip>
#include <iostream>
#include <thread>

#include "SolverConduction.hpp"

#define BCCAST static_cast<int16_t>

using Clock = std::chrono::high_resolution_clock;
using TimeType = std::chrono::time_point<std::chrono::high_resolution_clock>;

static void displayDT(TimeType startTime, TimeType endTime, std::string text)
{
    auto ellapsedTimeMeasure = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << text << static_cast<double>(ellapsedTimeMeasure.count())/1000.0 << " s" << std::endl;
}


SolverConduction::SolverConduction(const SolverConductionCreateInfo& solverCondInfos) :
Solver(solverCondInfos.solverInfos),
m_rho(solverCondInfos.rho),
m_k(solverCondInfos.k),
m_cv(solverCondInfos.cv)
{
    m_solverType = SOLVER_TYPE::Conduction;
    unsigned short dim = m_mesh.getDim();

    m_statesNumber = 1;

    m_mesh.setStatesNumber(1);

    m_lua["rho"] = m_rho;
    m_lua["k"] = m_k;
    m_lua["cv"] = m_cv;

    m_MPrev.resize((dim + 1), (dim + 1)); m_MPrev.setZero();
    for(unsigned short k = 0 ; k < m_N2.size() ; ++k)
    {
        m_MPrev += (m_N2[k].topLeftCorner(1, dim + 1)).transpose()*m_N2[k].topLeftCorner(1, dim + 1)*m_w2[k];
    }
    m_MPrev *= m_mesh.getRefElementSize(dim)*m_rho*m_cv;

    setInitialCondition();

    for(std::size_t n = 0 ; n < m_mesh.getNodesCount() ; ++n)
    {
        if(!m_mesh.getNode(n).isBound())
            continue;

        std::string bcName = m_mesh.getNodeType(n);

        bool boundaryT = false;
        bool boundaryQ = false;

        const auto& bcFunc1 = m_lua[bcName + "T"];
        boundaryT = bcFunc1.valid();

        const auto& bcFunc2 = m_lua[bcName + "Q"];
        boundaryQ = bcFunc2.valid();

        if(boundaryT && boundaryQ)
            throw std::runtime_error("Cannot apply a BC on T and Q on the same boundary !");
        else if(boundaryT)
            m_mesh.setUserDefTag(n, BCCAST(BC::T));
        else if(boundaryQ)
            m_mesh.setUserDefTag(n, BCCAST(BC::Q));
        else
            m_mesh.setUserDefTag(n, BCCAST(BC::None));
    }
}

void SolverConduction::applyBoundaryConditions(Eigen::VectorXd& b, const Eigen::VectorXd& qPrev)
{
    assert(m_mesh.getNodesCount() != 0);

    const uint8_t dim = m_mesh.getDim();
    const std::size_t nodesCount = m_mesh.getNodesCount();
    const std::size_t facetsCount = m_mesh.getFacetsCount();
    const std::size_t noPerFacet = m_mesh.getFacet(0).getNodesCount();

    for(std::size_t f = 0 ; f < facetsCount ; ++f)
    {
        const Facet& facet = m_mesh.getFacet(f);
        bool boundaryQ = true;
        for(uint8_t n = 0 ; n < noPerFacet ; ++n)
        {
            if(m_mesh.getNode(facet.getNodeIndex(n)).getUserDefTag() != BCCAST(BC::Q))
            {
                boundaryQ = false;
                break;
            }
        }
        if(!boundaryQ)
            continue;


        std::vector<double> t(noPerFacet);
        for(uint8_t n = 0 ; n < noPerFacet ; ++n)
        {
            const Node& node = m_mesh.getNode(facet.getNodeIndex(n));

            std::vector<double> result;
            result = m_lua[m_mesh.getNodeType(facet.getNodeIndex(n)) + "Q"](node.getPosition(),
                                                        m_mesh.getBoundNodeInitPos(n),
                                                        m_currentTime + m_currentDT).get<std::vector<double>>();

            t[n] = result[0];
        }

        //This is bad as it do not evalute q at gauss points
        auto evaluateTalongEdge = [&](double xi) -> double
        {
            double edgeLength = facet.getDetJ()*m_mesh.getRefElementSize(1);

            double x = facet.getJ(0, 0)*xi + edgeLength/2;

            return (t[1] - t[0])/edgeLength*x + t[0];
        };

        std::vector<std::array<double, 3>> gp = m_mesh.getGaussPoints(1, 2);

        //Ff = S Nv^T t dS
        Eigen::VectorXd Ff(dim*noPerFacet); Ff.setZero();
        for (unsigned short k = 0 ; k < m_ldN2.size() ; ++k)
        {
            Eigen::VectorXd Ffk (dim*noPerFacet);
            for(unsigned short i = 0 ; i < noPerFacet ; ++i)
            {
                Ffk[i] = m_ldN2[k](i)*evaluateTalongEdge(gp[k][0]);
            }

            Ff += Ffk*m_ldw2[k];
        }

        Ff *= m_mesh.getRefElementSize(dim - 1)*facet.getDetJ();

        for(unsigned short i = 0 ; i < noPerFacet ; ++i)
        {
            b(facet.getNodeIndex(i)) -= Ff[i]; //See Galerkin formulation
            //std::cout << Ff[i] << std::endl;
        }
    }

    //Do not parallelize this
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_mesh.getNode(n);
        if(node.getUserDefTag() == BCCAST(BC::T))
        {
            std::vector<double> result;
            result = m_lua[m_mesh.getNodeType(n) + "T"](node.getPosition(),
                                                        m_mesh.getBoundNodeInitPos(n),
                                                        m_currentTime + m_currentDT).get<std::vector<double>>();

            b(n) = result[0];
        }
    }
}

void SolverConduction::displaySolverParams() const noexcept
{
    std::cout << "Initial nodes number: " << m_mesh.getNodesCount() << std::endl;
    std::cout << "Initial elements number: " << m_mesh.getElementsCount() << std::endl;
    m_mesh.displayToConsole();
    std::cout << "Eigen sparse solver: SparseLU" << std::endl;
#if defined(_OPENMP)
    std::cout << "Number of OpenMP threads: " << m_numOMPThreads << "/" << omp_get_num_procs() << std::endl;
#endif
    std::cout << "Density: " << m_rho << " kg/m^3" << std::endl;
    std::cout << "Thermal conductivity: " << m_k << " W/(m K)" << std::endl;
    std::cout << "Specific heat capacity: " << m_cv << " J/(kg K)" << std::endl;
    std::cout << "End simulation time: " << m_endTime << " s" << std::endl;
    std::cout << "Adapt time step: " << (m_adaptDT ? "yes" : "no") << std::endl;
    if(m_adaptDT)
    {
        std::cout << "Maximum time step: " << m_maxDT << " s" << std::endl;
    }
    else
        std::cout << "Time step: " << m_maxDT << " s" << std::endl;
}

void SolverConduction::solveProblem(bool verboseOutput)
{
    std::cout   << "================================================================"
                << std::endl
                << "                EXECUTING PFEM CONDUCTION SOLVER                "
                << std::endl
                << "================================================================"
                << std::endl;

    displaySolverParams();

    std::cout << "----------------------------------------------------------------" << std::endl;

    for(auto& extractor : m_pExtractor)
    {
        extractor->update();
    }

    while(m_currentTime < m_endTime)
    {
        if(verboseOutput)
        {
            std::cout << "----------------------------------------------------------------" << std::endl;
            std::cout << "Solving time step: " << m_currentTime + m_currentDT
                      << "/" << m_endTime << " s, dt = " << m_currentDT << std::endl;
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
        else
        {
            std::cout << std::fixed << std::setprecision(3);
            std::cout << "\r" << "Solving time step: " << m_currentTime + m_currentDT
                      << "/" << m_endTime << " s, dt = ";
            std::cout << std::scientific;
            std::cout << m_currentDT << " s" << std::flush;
        }

        if(solveCurrentTimeStep(verboseOutput))
        {
            for(auto& extractor : m_pExtractor)
            {
                extractor->update();
            }
        }
        else
        {
            throw std::runtime_error("The solver does not seem to converge");
        }
    }

    std::cout << std::endl;
}

bool SolverConduction::solveCurrentTimeStep(bool verboseOutput)
{
    TimeType startTime, endTime;
    startTime = Clock::now();

    Eigen::VectorXd qPrev = getQFromNodesStates(0, m_statesNumber - 1);    //The precedent solution

    Eigen::SparseMatrix<double> A(m_mesh.getNodesCount(), m_mesh.getNodesCount());    //The matrix A representing the problem: [M/dt+K -D^T; C/dt-D L].
    Eigen::VectorXd b(m_mesh.getNodesCount());                                          //The vector b representing the problem: [M/dt*qPrev + F; H].

    Eigen::VectorXd qIter(m_mesh.getNodesCount()); qIter.setZero();

    buildPicardSystem(A, b, qPrev);

    applyBoundaryConditions(b, qPrev);

    m_solverLU.compute(A);
    if(m_solverLU.info() == Eigen::Success)
    {
        qIter = m_solverLU.solve(b);
    }
    else
    {
        m_mesh.restoreNodesList();
        return false;
    }

    setNodesStatesfromQ(qIter, 0, m_statesNumber - 1);

    m_currentTime += m_currentDT;
    m_currentStep++;

    endTime = Clock::now();
    if(verboseOutput)
        displayDT(startTime, endTime, "Picard algorithm solved in ");

    startTime = Clock::now();
    m_mesh.remesh(verboseOutput);
    endTime = Clock::now();
    if(verboseOutput)
        displayDT(startTime, endTime, "Remeshing done in ");

    return true;
}

void SolverConduction::buildPicardSystem(Eigen::SparseMatrix<double>& A,
                                         Eigen::VectorXd& b,
                                         const Eigen::VectorXd& qPrev)
{
    const unsigned short dim = m_mesh.getDim();
    const unsigned short noPerEl = dim + 1;
    const unsigned int tripletPerElm = (2*noPerEl*noPerEl);
    const std::size_t nElm = m_mesh.getElementsCount();
    const std::size_t nNodes = m_mesh.getNodesCount();

    std::vector<std::pair<std::size_t, double>> bDoubletList(noPerEl*nElm);
    std::vector<Eigen::Triplet<double>> indexA(tripletPerElm*nElm);

    b.setZero();

    Eigen::MatrixXd MPrev = m_MPrev/m_currentDT;

    Eigen::setNbThreads(1);
    //#pragma omp parallel for default(shared))
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_mesh.getElement((elm));

        Eigen::MatrixXd Be = getB(element);
        Eigen::MatrixXd Bep(dim, dim + 1);
        if(dim == 2)
            Bep << Be.block(0, 0, 1, 3), Be.block(1, 3, 1, 3);
        else
            Bep << Be.block(0, 0, 1, 4), Be.block(1, 4, 1, 4), Be.block(2, 8, 1, 4);

        //Me = S rho Nv^T Nv dV
        Eigen::MatrixXd Me = MPrev*element.getDetJ();

        Eigen::MatrixXd Le = m_mesh.getRefElementSize(dim)*m_k*Bep.transpose()*Bep*element.getDetJ();

        Eigen::VectorXd thetaPrev(noPerEl);
        if(dim == 2)
            thetaPrev << qPrev[element.getNodeIndex(0)], qPrev[element.getNodeIndex(1)], qPrev[element.getNodeIndex(2)];
        else
            thetaPrev << qPrev[element.getNodeIndex(0)], qPrev[element.getNodeIndex(1)], qPrev[element.getNodeIndex(2)], qPrev[element.getNodeIndex(3)];

        Eigen::VectorXd MThetae = Me * thetaPrev;

        std::size_t countA = 0;
        std::size_t countb = 0;

        for(unsigned short i = 0 ; i < (dim+1) ; ++i)
        {
            const Node& ni = m_mesh.getNode(element.getNodeIndex(i));

            for(unsigned short j = 0 ; j < (dim+1) ; ++j)
            {
                if(ni.getUserDefTag() != BCCAST(BC::T))
                {
                    indexA[tripletPerElm*elm + countA] =
                        Eigen::Triplet<double>(element.getNodeIndex(i),
                                               element.getNodeIndex(j),
                                               Me(i, j));

                    countA++;

                    indexA[tripletPerElm*elm + countA] =
                        Eigen::Triplet<double>(element.getNodeIndex(i),
                                               element.getNodeIndex(j),
                                               Le(i, j));

                    countA++;
                }
            }

            /************************************************************************
                                            Build h
            ************************************************************************/
            bDoubletList[noPerEl*elm + countb] = std::make_pair(element.getNodeIndex(i), MThetae[i]);

            countb++;
        }
    }
    Eigen::setNbThreads(m_numOMPThreads);

    //Best would be to know the number of nodes in which case :/
    //This can still be fasten using openmp but will never be as good as using []
    //with preallocated memory
    for(std::size_t n = 0 ; n < nNodes ; ++n)
    {
        const Node& node = m_mesh.getNode(n);

        if(node.getUserDefTag() == BCCAST(BC::T))
        {
            indexA.push_back(Eigen::Triplet<double>(n, n, 1));
        }
    }

    /********************************************************************************
                                        Compute A and b
    ********************************************************************************/
    A.setFromTriplets(indexA.begin(), indexA.end());

    for(auto& doublet : bDoubletList)
    {
        b[doublet.first] += doublet.second;
    }
}
