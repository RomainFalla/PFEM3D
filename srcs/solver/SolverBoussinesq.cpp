#include <chrono>
#include <iomanip>
#include <iostream>

#include "SolverBoussinesq.hpp"

#define BCCAST static_cast<int16_t>

using Clock = std::chrono::high_resolution_clock;
using TimeType = std::chrono::time_point<std::chrono::high_resolution_clock>;

static void displayDT(TimeType startTime, TimeType endTime, std::string text)
{
    auto ellapsedTimeMeasure = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << text << static_cast<double>(ellapsedTimeMeasure.count())/1000.0 << " s" << std::endl;
}


SolverBoussinesq::SolverBoussinesq(const SolverBoussinesqCreateInfo& solverBoussinesqInfos) :
Solver(solverBoussinesqInfos.solverInfos),
m_rho(solverBoussinesqInfos.rho),
m_mu(solverBoussinesqInfos.mu),
m_k(solverBoussinesqInfos.k),
m_cv(solverBoussinesqInfos.cv),
m_alpha(solverBoussinesqInfos.alpha),
m_Tr(solverBoussinesqInfos.Tr),
m_picardRelTol(solverBoussinesqInfos.picardRelTol),
m_picardMaxIter(solverBoussinesqInfos.picardMaxIter),
m_coeffDTincrease(solverBoussinesqInfos.coeffDTincrease),
m_coeffDTdecrease(solverBoussinesqInfos.coeffDTdecrease)
{
    m_solverType = SOLVER_TYPE::Boussinesq;
    unsigned short dim = m_mesh.getDim();

    m_statesNumber = dim + 2; //(u,v, (w, ) p, T)

    m_mesh.setStatesNumber(m_statesNumber);

    m_lua["rho"] = m_rho;
    m_lua["mu"] = m_mu;
    m_lua["k"] = m_k;
    m_lua["cv"] = m_cv;
    m_lua["alpha"] = m_alpha;
    m_lua["Tr"] = m_Tr;

    m_MPrev.resize(dim*(dim + 1), dim*(dim + 1)); m_MPrev.setZero();
    for(unsigned short k = 0 ; k < m_N2.size() ; ++k)
    {
        m_MPrev += m_N2[k].transpose()*m_N2[k]*m_w2[k];
    }
    m_MPrev *= m_mesh.getRefElementSize(dim)*m_rho;

    m_MPrevT.resize((dim + 1), (dim + 1)); m_MPrevT.setZero();
    for(unsigned short k = 0 ; k < m_N2.size() ; ++k)
    {
        m_MPrevT += (m_N2[k].topLeftCorner(1, dim + 1)).transpose()*m_N2[k].topLeftCorner(1, dim + 1)*m_w2[k];
    }
    m_MPrevT *= m_mesh.getRefElementSize(dim)*m_rho*m_cv;

    if(dim == 2)
    {
        m_ddev.resize(3, 3);
        m_ddev << 2, 0, 0,
                  0, 2, 0,
                  0, 0, 1;
    }
    else
    {
        m_ddev.resize(6, 6);
        m_ddev << 2, 0, 0, 0, 0, 0,
                  0, 2, 0, 0, 0, 0,
                  0, 0, 2, 0, 0, 0,
                  0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 1;
    }

    m_ddev *= m_mu;

    setInitialCondition();

    for(std::size_t n = 0 ; n < m_mesh.getNodesCount() ; ++n)
    {
        if(!m_mesh.getNode(n).isBound())
            continue;

        std::string bcName = m_mesh.getNodeType(n);

        bool boundaryT = false;
        bool boundaryQ = false;
        bool boundaryV = false;

        const auto& bcFunc1 = m_lua[bcName + "T"];
        boundaryT = bcFunc1.valid();

        const auto& bcFunc2 = m_lua[bcName + "Q"];
        boundaryQ = bcFunc2.valid();

        const auto& bcFunc3 = m_lua[bcName + "V"];
        boundaryV = bcFunc3.valid();

        if(boundaryT && boundaryQ)
            throw std::runtime_error("Cannot apply a BC on T and Q on the same boundary !");
        else if(boundaryT && boundaryV)
            m_mesh.setUserDefTag(n, BCCAST(BC::TandV));
        else if(boundaryQ && boundaryV)
            m_mesh.setUserDefTag(n, BCCAST(BC::QandV));
        else if(boundaryT)
            m_mesh.setUserDefTag(n, BCCAST(BC::T));
        else if(boundaryQ)
            m_mesh.setUserDefTag(n, BCCAST(BC::Q));
        else if(boundaryV)
            m_mesh.setUserDefTag(n, BCCAST(BC::V));
        else
            m_mesh.setUserDefTag(n, BCCAST(BC::None));
    }

    if(m_gamma < 1e-15)
        m_mesh.setComputeNormalCurvature(false);
}

void SolverBoussinesq::applyBoundaryConditionsHeat(Eigen::VectorXd& b, const Eigen::VectorXd& qPrev)
{
    assert(m_mesh.getNodesCount() != 0);

    const uint8_t dim = m_mesh.getDim();
    const std::size_t nodesCount = m_mesh.getNodesCount();
    const std::size_t facetsCount = m_mesh.getFacetsCount();
    const std::size_t noPerFacet = m_mesh.getFacet(0).getNodesCount();

//    for(std::size_t f = 0 ; f < facetsCount ; ++f)
//    {
//        const Facet& facet = m_mesh.getFacet(f);
//        bool boundaryQ = true;
//        for(uint8_t n = 0 ; n < noPerFacet ; ++n)
//        {
//            if(m_mesh.getNode(facet.getNodeIndex(n)).getUserDefTag() != BCCAST(BC::Q) &&
//               m_mesh.getNode(facet.getNodeIndex(n)).getUserDefTag() != BCCAST(BC::QandV))
//            {
//                boundaryQ = false;
//                break;
//            }
//        }
//        if(!boundaryQ)
//            continue;
//
//
//        std::vector<double> t(noPerFacet);
//        for(uint8_t n = 0 ; n < noPerFacet ; ++n)
//        {
//            const Node& node = m_mesh.getNode(facet.getNodeIndex(n));
//
//            std::vector<double> result;
//            result = m_lua[m_mesh.getNodeType(facet.getNodeIndex(n)) + "Q"](node.getPosition(),
//                                                        m_mesh.getBoundNodeInitPos(n),
//                                                        m_currentTime + m_currentDT).get<std::vector<double>>();
//
//            t[n] = result[0];
//        }
//
//        //This is bad as it do not evalute q at gauss points
//        auto evaluateTalongEdge = [&](double xi) -> double
//        {
//            double edgeLength = facet.getDetJ()*m_mesh.getRefElementSize(1);
//
//            double x = facet.getJ(0, 0)*xi + edgeLength/2;
//
//            return (t[1] - t[0])/edgeLength*x + t[0];
//        };
//
//        std::vector<std::array<double, 3>> gp = m_mesh.getGaussPoints(1, 2);
//
//        //Ff = S Nv^T t dS
//        Eigen::VectorXd Ff(dim*noPerFacet); Ff.setZero();
//        for (unsigned short k = 0 ; k < m_ldN2.size() ; ++k)
//        {
//            Eigen::VectorXd Ffk (dim*noPerFacet);
//            for(unsigned short i = 0 ; i < noPerFacet ; ++i)
//            {
//                Ffk[i] = m_ldN2[k](i)*evaluateTalongEdge(gp[k][0]);
//            }
//
//            Ff += Ffk*m_ldw2[k];
//        }
//
//        Ff *= m_mesh.getRefElementSize(dim - 1)*facet.getDetJ();
//
//        for(unsigned short i = 0 ; i < noPerFacet ; ++i)
//        {
//            b(facet.getNodeIndex(i)) -= Ff[i]; //See Galerkin formulation
//            //std::cout << Ff[i] << std::endl;
//        }
//    }

    //Do not parallelize this
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_mesh.getNode(n);
        if(node.isFree())
        {
            b(n) = qPrev[n];
        }
        else if(node.getUserDefTag() == BCCAST(BC::T) ||
                node.getUserDefTag() == BCCAST(BC::TandV))
        {
            std::vector<double> result;
            result = m_lua[m_mesh.getNodeType(n) + "T"](node.getPosition(),
                                                        m_mesh.getBoundNodeInitPos(n),
                                                        m_currentTime + m_currentDT).get<std::vector<double>>();

            b(n) = result[0];
        }
    }
}

void SolverBoussinesq::applyBoundaryConditionsMom(Eigen::VectorXd& b, const Eigen::VectorXd& qPrev)
{
    assert(m_mesh.getNodesCount() != 0);
    const uint8_t dim = m_mesh.getDim();
    const std::size_t nodesCount = m_mesh.getNodesCount();

    const std::size_t facetsCount = m_mesh.getFacetsCount();
    const std::size_t noPerFacet = m_mesh.getFacet(0).getNodesCount();

    for(std::size_t f = 0 ; f < facetsCount ; ++f)
    {
        if(m_gamma < 1e-15)
            continue;

        const Facet& facet = m_mesh.getFacet(f);

        bool onFS = true;
        for(unsigned short n = 0 ; n < noPerFacet ; ++n)
        {
            if(!m_mesh.getNode(facet.getNodeIndex(n)).isOnFreeSurface())
            {
                onFS = false;
                break;
            }
        }
        if(!onFS)
            continue;

        if(dim == 3)
        {
            std::cerr << "Surface tension boundary condition not implemented in 3D" << std::endl;
            break;
        }

        std::vector<Eigen::VectorXd> t(noPerFacet);
        for(uint8_t n = 0 ; n < noPerFacet ; ++n)
        {
            double curvature = m_mesh.getFreeSurfaceCurvature(facet.getNodeIndex(n));
            std::array<double, 3> normal = m_mesh.getFreeSurfaceNormal(facet.getNodeIndex(n));


            t[n].resize(dim);
            for(uint8_t d = 0 ; d < dim ; ++d)
                t[n](d) = m_gamma*curvature*normal[d];
        }

        auto evaluateTalongEdge = [&](double xi, uint8_t coordinate) -> double
        {
            double edgeLength = facet.getDetJ()*m_mesh.getRefElementSize(1);

            double x = facet.getJ(0, 0)*xi + edgeLength/2;

            return (t[1][coordinate] - t[0][coordinate])/edgeLength*x + t[0][coordinate];

        };

        std::vector<std::array<double, 3>> gp = m_mesh.getGaussPoints(1, 2);

        //Ff = S Nv^T t dS
        Eigen::VectorXd Ff(dim*noPerFacet); Ff.setZero();
        for (unsigned short k = 0 ; k < m_ldN2.size() ; ++k)
        {
            Eigen::VectorXd Ffk (dim*noPerFacet);
            for(unsigned short i = 0 ; i < noPerFacet ; ++i)
            {
                for(unsigned short d = 0 ; d < dim ; ++d)
                {
                    Ffk[i + d*noPerFacet] = m_ldN2[k](i)*evaluateTalongEdge(gp[k][0], d);
                }
            }

            Ff += Ffk*m_ldw2[k];
        }

        Ff *= m_mesh.getRefElementSize(dim - 1)*facet.getDetJ();

        for(unsigned short i = 0 ; i < noPerFacet ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                //std::cout << Ff[d*noPerFace + i] << std::endl;
                b(facet.getNodeIndex(i) + d*nodesCount) += Ff[d*noPerFacet + i];
            }
        }
    }

    //Do not parallelize this
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_mesh.getNode(n);
        if(node.isFree())
        {
            b(n + dim*nodesCount) = 0;

            if(!node.isBound())
            {
                for(uint8_t d = 0 ; d < dim - 1 ; ++d)
                {
                    b(n + d*nodesCount) = qPrev(n + d*nodesCount);
                }

                b(n + (dim - 1)*nodesCount) = qPrev(n + (dim - 1)*nodesCount) - m_currentDT*m_gravity;
            }
        }

        if(node.getUserDefTag() == BCCAST(BC::V) ||
           node.getUserDefTag() == BCCAST(BC::QandV) ||
           node.getUserDefTag() == BCCAST(BC::TandV))
        {
            std::vector<double> result;
            result = m_lua[m_mesh.getNodeType(n) + "V"](node.getPosition(),
                                                  m_mesh.getBoundNodeInitPos(n),
                                                  m_currentTime + m_currentDT).get<std::vector<double>>();

            for(uint8_t d = 0 ; d < dim ; ++d)
            {
                b(n + d*nodesCount) = result[d];
            }
        }

        if(m_strongPAtFS && node.isOnFreeSurface())
        {
            b(n + dim*nodesCount) = 0;
        }
    }
}

void SolverBoussinesq::displaySolverParams() const noexcept
{
    std::cout << "Initial nodes number: " << m_mesh.getNodesCount() << "\n";
    std::cout << "Initial elements number: " << m_mesh.getElementsCount() << "\n";
    m_mesh.displayToConsole();
    std::cout << "Eigen sparse solver: SparseLU" << "\n";
#if defined(_OPENMP)
    std::cout << "Number of OpenMP threads: " << m_numOMPThreads << "/" << omp_get_num_procs() << "\n";
#endif
    std::cout << "Gravity: " << m_gravity << " m/s^2" << "\n";
    std::cout << "Density: " << m_rho << " kg/m^3" << "\n";
    std::cout << "Viscosity: " << m_mu << " Pa s" << "\n";
    std::cout << "Surface tension: " << m_gamma << " N/m" << "\n";
    std::cout << "Thermal conductivity: " << m_k << " W/(m K)" << "\n";
    std::cout << "Specific heat capacity: " << m_cv << " W/(m K)" << "\n";
    std::cout << "Picard relative tolerance: " << m_picardRelTol  << "\n";
    std::cout << "Maximum picard iteration number: " << m_picardMaxIter << "\n";
    std::cout << "End simulation time: " << m_endTime << " s" << "\n";
    std::cout << "Adapt time step: " << (m_adaptDT ? "yes" : "no") << "\n";
    if(m_adaptDT)
    {
        std::cout << "Maximum time step: " << m_maxDT << " s" << "\n";
        std::cout << "Time step reduction coefficient: " << m_coeffDTdecrease << "\n";
        std::cout << "Time step increase coefficient: " << m_coeffDTincrease << "\n";
    }
    else
        std::cout << "Time step: " << m_maxDT << " s" << "\n";

    std::cout << std::flush;
}

void SolverBoussinesq::solveProblem(bool verboseOutput)
{
    std::cout   << "================================================================"
                << "\n"
                << "                 EXECUTING PFEM BOUSSINESQ SOLVER               "
                << "\n"
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

            if(m_adaptDT)
            {
                if(m_picardCurrentNumIter < m_picardMaxIter)
                    m_currentDT = std::min(m_maxDT, m_currentDT*m_coeffDTincrease);

            }
        }
        else
        {
            if(m_adaptDT)
            {
                if(verboseOutput)
                {
                    std::cout << "* The algorithm did not converge - Reducing the time step"
                              << std::endl;
                }

                m_currentDT /= m_coeffDTdecrease;
            }
            else
            {
                throw std::runtime_error("The solver does not seem to converge");
            }
        }
    }

    std::cout << std::endl;
}

bool SolverBoussinesq::solveCurrentTimeStep(bool verboseOutput)
{
    TimeType startTime, endTime, startTimeMeasure, endTimeMeasure;
    startTime = Clock::now();

    const unsigned short dim = m_mesh.getDim();

    /***************************
        Solve the heat equation
    ****************************/

    if(verboseOutput)
        std::cout << "Solving heat equation" << std::endl;

    Eigen::VectorXd qPrevHeat = getQFromNodesStates(m_statesNumber - 1, m_statesNumber - 1);

    Eigen::SparseMatrix<double> Aheat(m_mesh.getNodesCount(), m_mesh.getNodesCount());
    Eigen::VectorXd bHeat(m_mesh.getNodesCount());

    m_mesh.saveNodesList();

    Eigen::VectorXd qHeat(m_mesh.getNodesCount()); qHeat.setZero();

    buildAHeat(Aheat, bHeat, qPrevHeat);
    applyBoundaryConditionsHeat(bHeat, qPrevHeat);

    m_solverLU.compute(Aheat);
    if(m_solverLU.info() == Eigen::Success)
    {
        qHeat = m_solverLU.solve(bHeat);
    }
    else
    {
        m_mesh.restoreNodesList();
        return false;
    }

    setNodesStatesfromQ(qHeat, m_statesNumber - 1, m_statesNumber - 1);

    /********************************
        Solve the momentum equation
    *********************************/

    if(verboseOutput)
        std::cout << "Solving momentum equation" << std::endl;

    Eigen::VectorXd qPrevMom = getQFromNodesStates(0, m_statesNumber - 2);

    Eigen::SparseMatrix<double> AMom((m_statesNumber-1)*m_mesh.getNodesCount(), (m_statesNumber-1)*m_mesh.getNodesCount());
    Eigen::VectorXd bMom((m_statesNumber-1)*m_mesh.getNodesCount());

    m_mesh.saveNodesList();

    m_picardCurrentNumIter = 0;
    Eigen::VectorXd qIterMom((m_statesNumber-1)*m_mesh.getNodesCount()); qIterMom.setZero();
    Eigen::VectorXd qIterPrevMom((m_statesNumber-1)*m_mesh.getNodesCount()); qIterPrevMom.setZero();
    double res = std::numeric_limits<double>::max();

    while(res > m_picardRelTol)
    {
        if(verboseOutput)
        {
            std::cout << "Picard algorithm (mesh position) - iteration ("
                      << m_picardCurrentNumIter << ")" << std::endl;
        }

        qIterPrevMom = qIterMom;

        buildPicardSystemMom(AMom, bMom, qPrevMom);
        applyBoundaryConditionsMom(bMom, qPrevMom);

        m_solverLU.compute(AMom);
        if(m_solverLU.info() == Eigen::Success)
        {
            qIterMom= m_solverLU.solve(bMom);
        }
        else
        {
            m_mesh.restoreNodesList();
            return false;
        }

        setNodesStatesfromQ(qIterMom, 0, m_statesNumber - 2);
        Eigen::VectorXd deltaPos = qIterMom*m_currentDT;
        m_mesh.updateNodesPositionFromSave(std::vector<double> (deltaPos.data(), deltaPos.data() + dim*m_mesh.getNodesCount()));
        m_picardCurrentNumIter++;

        if(m_picardCurrentNumIter == 1)
            res = std::numeric_limits<double>::max();
        else
        {
            double num{0};
            double den{0};

            for(std::size_t n = 0 ; n < m_mesh.getNodesCount() ; ++n)
            {
                const Node& node = m_mesh.getNode(n);

                if(!node.isFree())
                {
                    for(unsigned short d = 0 ; d < dim ; ++d)
                    {
                        num += (qIterMom(n + d*m_mesh.getNodesCount()) - qIterPrevMom(n + d*m_mesh.getNodesCount()))*(qIterMom(n + d*m_mesh.getNodesCount()) - qIterPrevMom(n + d*m_mesh.getNodesCount()));
                        den += qIterPrevMom(n + d*m_mesh.getNodesCount())*qIterPrevMom(n + d*m_mesh.getNodesCount());
                    }
                }
            }

            res = std::sqrt(num/den);
        }

        if(verboseOutput)
        {
            startTimeMeasure = Clock::now();
            std::cout << " * Relative 2-norm of q: " << res << " vs "
                      << m_picardRelTol << std::endl;
        }

        if(m_picardCurrentNumIter >= m_picardMaxIter || std::isnan(res))
        {
            m_mesh.restoreNodesList();
            return false;
        }
    }

    m_currentTime += m_currentDT;
    m_currentStep++;

    endTime = Clock::now();
    if(verboseOutput)
        displayDT(startTime, endTime, "Problem solved in ");

    startTime = Clock::now();
    m_mesh.remesh(verboseOutput);
    endTime = Clock::now();
    if(verboseOutput)
        displayDT(startTime, endTime, "Remeshing done in ");

    return true;
}

void SolverBoussinesq::buildAHeat(Eigen::SparseMatrix<double>& A,
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

    Eigen::MatrixXd MPrev = m_MPrevT/m_currentDT;

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
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
                if(ni.getUserDefTag() != BCCAST(BC::T) && ni.getUserDefTag() != BCCAST(BC::TandV))
                {
                    if(!ni.isFree())
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(element.getNodeIndex(i),
                                                   element.getNodeIndex(j),
                                                   Me(i, j));
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

        if(node.getUserDefTag() == BCCAST(BC::T) || node.getUserDefTag() == BCCAST(BC::TandV) || node.isFree())
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

void SolverBoussinesq::buildPicardSystemMom(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b,
                                             const Eigen::VectorXd& qPrev)
{
    const uint8_t dim = m_mesh.getDim();
    const uint8_t noPerEl = m_mesh.getElement(0).getNodeCount();
    const uint16_t tripletPerElm = (dim*noPerEl*noPerEl + dim*noPerEl*dim*noPerEl + 3*noPerEl*dim*noPerEl + noPerEl*noPerEl);
    const uint16_t doubletPerElm = (2*dim*noPerEl + 2*noPerEl);
    const std::size_t nElm = m_mesh.getElementsCount();
    const std::size_t nNodes = m_mesh.getNodesCount();

    std::vector<Eigen::Triplet<double>> indexA(tripletPerElm*nElm);
    std::vector<std::pair<std::size_t, double>> indexb(doubletPerElm*nElm); b.setZero();

    Eigen::MatrixXd MPrev = m_MPrev/m_currentDT;

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_mesh.getElement(elm);

        Eigen::MatrixXd Be = getB(element);
        Eigen::MatrixXd Bep(dim, dim + 1);
        if(dim == 2)
            Bep << Be.block(0, 0, 1, 3), Be.block(1, 3, 1, 3);
        else
            Bep << Be.block(0, 0, 1, 4), Be.block(1, 4, 1, 4), Be.block(2, 8, 1, 4);

        double tauPSPG = computeTauPSPG(element);

        //Me = S rho Nv^T Nv dV
        Eigen::MatrixXd Me = MPrev*element.getDetJ();

        //Ke = S Bv^T ddev Bv dV
        Eigen::MatrixXd Ke = m_mesh.getRefElementSize(dim)*Be.transpose()*m_ddev*Be*element.getDetJ(); //same matrices for all gauss point ^^

        //De = S Np^T m Bv dV
        Eigen::MatrixXd De((dim+1),dim*(dim+1)); De.setZero();

        for (unsigned short k = 0 ; k < m_N1.size() ; ++k)
            De += (m_N1[k].topLeftCorner(1, dim+1)).transpose()*m_m.transpose()*Be*m_w1[k];

        De *= m_mesh.getRefElementSize(dim)*element.getDetJ();

        //Ce = tauPSPG S Bp^T Nv dV
        Eigen::MatrixXd Ce((dim+1),dim*(dim+1)); Ce.setZero();

        for (unsigned short k = 0 ; k < m_N1.size() ; ++k)
            Ce += Bep.transpose()*m_N1[k]*m_w1[k];

        Ce *= (m_mesh.getRefElementSize(dim)*tauPSPG/m_currentDT)*element.getDetJ();

        //Le = tauPSPG S Bp^T Bp dV
        Eigen::MatrixXd Le = m_mesh.getRefElementSize(dim)*(tauPSPG/m_rho)*
                             Bep.transpose()*Bep*element.getDetJ();

        Eigen::VectorXd Telm = getElementState(element, m_statesNumber - 1);

        //Fe = S Nv^T bodyforce dV
        Eigen::VectorXd Fe(dim*noPerEl); Fe.setZero();
        for (unsigned short k = 0 ; k < m_N2.size() ; ++k)
        {
            double T = (m_N2[k].topLeftCorner(1, dim + 1)*Telm).value();
            Fe += (1-m_alpha*(T - m_Tr))*m_N2[k].transpose()*m_bodyForces*m_w2[k];
        }

        //Oh could I have forgot this ?
        Fe *= m_rho*m_mesh.getRefElementSize(dim)*element.getDetJ();

        //He = tauPSPG S Bp^T bodyforce dV
        Eigen::VectorXd He(noPerEl); He.setZero();
        for (unsigned short k = 0 ; k < m_N1.size() ; ++k)
        {
            double T = (m_N1[k].topLeftCorner(1, dim + 1)*Telm).value();
            He += (1-m_alpha*(T - m_Tr))*Bep.transpose()*m_bodyForces*m_w1[k];
        }

        He *= m_mesh.getRefElementSize(dim)*tauPSPG*element.getDetJ();

        Eigen::VectorXd vPrev(noPerEl*dim);
        if(dim == 2)
            vPrev << qPrev[element.getNodeIndex(0)], qPrev[element.getNodeIndex(1)], qPrev[element.getNodeIndex(2)],
                     qPrev[element.getNodeIndex(0) + nNodes], qPrev[element.getNodeIndex(1) + nNodes], qPrev[element.getNodeIndex(2) + nNodes];
        else
            vPrev << qPrev[element.getNodeIndex(0)], qPrev[element.getNodeIndex(1)], qPrev[element.getNodeIndex(2)], qPrev[element.getNodeIndex(3)],
                     qPrev[element.getNodeIndex(0) + nNodes], qPrev[element.getNodeIndex(1) + nNodes], qPrev[element.getNodeIndex(2) + nNodes], qPrev[element.getNodeIndex(3) + nNodes],
                     qPrev[element.getNodeIndex(0) + 2*nNodes], qPrev[element.getNodeIndex(1) + 2*nNodes], qPrev[element.getNodeIndex(2) + 2*nNodes], qPrev[element.getNodeIndex(3) + 2*nNodes];

        Eigen::VectorXd MvPreve = Me*vPrev;
        Eigen::VectorXd CvPreve = Ce*vPrev;

        std::size_t countA = 0;
        std::size_t countb = 0;

        for(unsigned short i = 0 ; i < noPerEl ; ++i)
        {
            const Node& ni = m_mesh.getNode(element.getNodeIndex(i));

            for(unsigned short j = 0 ; j < noPerEl ; ++j)
            {
                const Node& nj = m_mesh.getNode(element.getNodeIndex(j));

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
                                                   Me(i + d*noPerEl, j + d*noPerEl));
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
                    if(!(ni.isFree() || (m_strongPAtFS && ni.isOnFreeSurface())))
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
                    if(!(ni.isFree() || (m_strongPAtFS && ni.isOnFreeSurface())))
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(element.getNodeIndex(i) + dim*nNodes,
                                                   element.getNodeIndex(j) + d*nNodes,
                                                   Ce(i, j + d*noPerEl));
                    }
                    countA++;
                }

                /********************************************************************
                                            Build L
                ********************************************************************/
                if(!(ni.isFree() || (m_strongPAtFS && ni.isOnFreeSurface())))
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

                indexb[doubletPerElm*elm + countb] = std::make_pair(element.getNodeIndex(i) + d*nNodes, MvPreve(i + d*noPerEl));
                countb++;
            }

            /************************************************************************
                                            Build h
            ************************************************************************/

            indexb[doubletPerElm*elm + countb] = std::make_pair(element.getNodeIndex(i) + dim*nNodes, He(i));
            countb++;

            indexb[doubletPerElm*elm + countb] = std::make_pair(element.getNodeIndex(i) + dim*nNodes, CvPreve(i));
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

        if(node.isFree() || (m_strongPAtFS && node.isOnFreeSurface()))
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

double SolverBoussinesq::computeTauPSPG(const Element& element)
{
    const double h = std::sqrt(m_mesh.getRefElementSize(m_mesh.getDim())*element.getDetJ()/M_PI);

    double U = 0;
    for (unsigned short n = 0 ; n < element.getNodeCount() ; ++n)
    {
        const Node& node = m_mesh.getNode(element.getNodeIndex(n));

        double nodeU = 0;
        for (unsigned short d = 0 ; d < m_mesh.getDim() ; ++d)
        {
            nodeU += node.getState(d)*node.getState(d);
        }
        U += std::sqrt(nodeU);
    }
    U /= (m_mesh.getDim() + 1);

    return 1/std::sqrt((2/m_currentDT)*(2/m_currentDT) + (2*U/h)*(2*U/h)
                        + 9*(4*m_mu/(h*h*m_rho))*(4*m_mu/(h*h*m_rho)));
}
