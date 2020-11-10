#include <chrono>
#include <iomanip>
#include <iostream>

#include "SolverIncompressible.hpp"

using Clock = std::chrono::high_resolution_clock;
using TimeType = std::chrono::time_point<std::chrono::high_resolution_clock>;

static void displayDT(TimeType startTime, TimeType endTime, std::string text)
{
    auto ellapsedTimeMeasure = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << text << static_cast<double>(ellapsedTimeMeasure.count())/1000.0 << " s" << std::endl;
}


SolverIncompressible::SolverIncompressible(const SolverIncompCreateInfo& solverIncompInfos) :
Solver(solverIncompInfos.solverInfos),
m_rho(solverIncompInfos.rho),
m_mu(solverIncompInfos.mu),
m_gamma(solverIncompInfos.gamma),
m_picardRelTol(solverIncompInfos.picardRelTol),
m_picardMaxIter(solverIncompInfos.picardMaxIter),
m_coeffDTincrease(solverIncompInfos.coeffDTincrease),
m_coeffDTdecrease(solverIncompInfos.coeffDTdecrease)
{
    m_solverType = SOLVER_TYPE::Incompressible_PSPG;
    unsigned short dim = m_mesh.getDim();

    m_statesNumber = dim + 1;

    m_mesh.setStatesNumber(m_statesNumber);

    m_lua["rho"] = m_rho;
    m_lua["mu"] = m_mu;
    m_lua["gamma"] = m_gamma;

    m_MPrev.resize(dim*(dim + 1), dim*(dim + 1)); m_MPrev.setZero();
    for(unsigned short k = 0 ; k < m_N2.size() ; ++k)
    {
        m_MPrev += m_N2[k].transpose()*m_N2[k]*m_w2[k];
    }
    m_MPrev *= m_mesh.getRefElementSize(dim)*m_rho;

    Eigen::MatrixXd usefulBloc(dim, dim); usefulBloc.setZero();
    Eigen::MatrixXd nullBloc(dim, dim); nullBloc.setZero();
    m_MPrevGamma.resize(dim*dim, dim*dim); m_MPrevGamma.setZero();
    for(unsigned short k = 0 ; k < m_ldN2.size() ; ++k)
    {
        usefulBloc += m_ldN2[k].transpose()*m_ldN2[k]*m_ldw2[k];
    }
    usefulBloc *= m_mesh.getRefElementSize(dim - 1)*m_gamma;
    if(dim == 2)
        m_MPrevGamma << usefulBloc, nullBloc, nullBloc, usefulBloc;
    else
    {
        m_MPrevGamma << usefulBloc, nullBloc, nullBloc,
                        nullBloc, usefulBloc, nullBloc,
                        nullBloc, nullBloc, usefulBloc;
    }

    std::cout << m_MPrevGamma << std::endl;

    //For now body forces constant so no problem
    m_FPrev.resize(dim*(dim + 1)); m_FPrev.setZero();
    for(unsigned short k = 0 ; k < m_N1.size() ; ++k)
    {
        m_FPrev += m_N1[k].transpose()*m_bodyForces*m_w1[k];
    }
    m_FPrev *= m_mesh.getRefElementSize(dim)*m_rho;

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

    if(m_gamma == 0)
        m_mesh.setComputeNormalCurvature(false);
}

void SolverIncompressible::applyBoundaryConditions(Eigen::VectorXd& b, const Eigen::VectorXd& qPrev)
{
    assert(m_mesh.getNodesCount() != 0);
    const uint8_t dim = m_mesh.getDim();
    const std::size_t nodesCount = m_mesh.getNodesCount();

    const std::size_t facetsCount = m_mesh.getFacetsCount();
    const std::size_t noPerFacet = m_mesh.getFacet(0).getNodesCount();

    //TO DO: paralelize this (no call to lua)
    for(std::size_t f = 0 ; f < facetsCount ; ++f)
    {
        if(m_gamma == 0)
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

        Eigen::VectorXd kappa_n(dim*noPerFacet);
        for(uint8_t n = 0 ; n < noPerFacet; ++n)
        {
            double curvature = m_mesh.getFreeSurfaceCurvature(facet.getNodeIndex(n));
            std::array<double, 3> normal = m_mesh.getFreeSurfaceNormal(facet.getNodeIndex(n));

            for(uint8_t d = 0 ; d < dim ; ++d)
                kappa_n(n + d*noPerFacet) = curvature*normal[d];
        }

        Eigen::VectorXd Ff = m_MPrevGamma*kappa_n*facet.getDetJ();

        for(unsigned short i = 0 ; i < noPerFacet ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
                b(facet.getNodeIndex(i) + d*nodesCount) += Ff[d*noPerFacet + i];
        }
    }

    //Do not parallelize this (lua)
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

        if(node.isBound())
        {
            std::vector<double> result;
            result = m_lua[m_mesh.getNodeType(n)](node.getPosition(),
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

void SolverIncompressible::displaySolverParams() const noexcept
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

void SolverIncompressible::solveProblem(bool verboseOutput)
{
    std::cout   << "================================================================"
                << "\n"
                << "            EXECUTING PFEM INCOMPRESSIBLE PSPG SOLVER           "
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

bool SolverIncompressible::solveCurrentTimeStep(bool verboseOutput)
{
    TimeType startTime, endTime, startTimeMeasure, endTimeMeasure;
    startTime = Clock::now();

    const unsigned short dim = m_mesh.getDim();

    Eigen::VectorXd qPrev = getQFromNodesStates(0, m_statesNumber - 1);

    Eigen::SparseMatrix<double, Eigen::RowMajor> A(m_statesNumber*m_mesh.getNodesCount(), m_statesNumber*m_mesh.getNodesCount());
    Eigen::VectorXd b(m_statesNumber*m_mesh.getNodesCount());

    m_mesh.saveNodesList();

    m_picardCurrentNumIter = 0;
    Eigen::VectorXd qIter(m_statesNumber*m_mesh.getNodesCount()); qIter.setZero();
    Eigen::VectorXd qIterPrev(m_statesNumber*m_mesh.getNodesCount()); qIterPrev.setZero();
    double res = std::numeric_limits<double>::max();

    while(res > m_picardRelTol)
    {
        if(verboseOutput)
        {
            std::cout << "Picard algorithm (mesh position) - iteration ("
                      << m_picardCurrentNumIter << ")" << std::endl;
        }

        qIterPrev = qIter;

        buildPicardSystem(A, b, qPrev);
        applyBoundaryConditions(b, qPrev);

        m_solverIt.compute(A);

        if(m_solverIt.info() == Eigen::Success)
        {
            qIter = m_solverIt.solveWithGuess(b, qIterPrev);
        }
        else
        {
            m_mesh.restoreNodesList();
            return false;
        }

        setNodesStatesfromQ(qIter, 0, m_statesNumber - 1);
        Eigen::VectorXd deltaPos = qIter*m_currentDT;
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
                        num += (qIter(n + d*m_mesh.getNodesCount()) - qIterPrev(n + d*m_mesh.getNodesCount()))*(qIter(n + d*m_mesh.getNodesCount()) - qIterPrev(n + d*m_mesh.getNodesCount()));
                        den += qIterPrev(n + d*m_mesh.getNodesCount())*qIterPrev(n + d*m_mesh.getNodesCount());
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
        displayDT(startTime, endTime, "Picard algorithm solved in ");

    startTime = Clock::now();
    m_mesh.remesh(verboseOutput);
    endTime = Clock::now();
    if(verboseOutput)
        displayDT(startTime, endTime, "Remeshing done in ");

    return true;
}

void SolverIncompressible::buildPicardSystem(Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& b,
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

        //Fe = S Nv^T bodyforce dV
        Eigen::MatrixXd Fe = m_FPrev*element.getDetJ();

        //He = tauPSPG S Bp^T bodyforce dV
        Eigen::MatrixXd He = m_mesh.getRefElementSize(dim)*tauPSPG*Bep.transpose()*
                             m_bodyForces*element.getDetJ();

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

double SolverIncompressible::computeTauPSPG(const Element& element)
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
