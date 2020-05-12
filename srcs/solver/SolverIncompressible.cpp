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


SolverIncompressible::SolverIncompressible(const nlohmann::json& j, const std::string& mshName) :
Solver(j, mshName)
{
    m_solverType = SOLVER_TYPE::Incompressible_PSPG;
    unsigned short dim = m_mesh.getDim();

    m_statesNumber = dim + 1;

    m_mesh.setStatesNumber(m_statesNumber);

    m_rho               = j["Solver"]["Fluid"]["rho"].get<double>();
    m_mu                = j["Solver"]["Fluid"]["mu"].get<double>();

    m_picardRelTol           = j["Solver"]["Picard"]["relTol"].get<double>();
    m_picardMaxIter          = j["Solver"]["Picard"]["maxIter"].get<unsigned int>();

    m_coeffDTincrease    = j["Solver"]["Time"]["coeffDTincrease"].get<double>();
    m_coeffDTdecrease    = j["Solver"]["Time"]["coeffDTdecrease"].get<double>();

    std::vector<std::string> whatCanBeWritten;
    if(dim == 2)
        whatCanBeWritten = {"u", "v", "p", "ke", "velocity"};
    else
        whatCanBeWritten = {"u", "v", "w", "p", "ke", "velocity"};

    auto extractors = j["Solver"]["Extractors"];
    unsigned short GMSHExtractorCount = 0;
    for(auto& extractor : extractors)
    {
        if(extractor["type"].get<std::string>() == "Point")
        {
            m_pExtractor.push_back(std::make_unique<PointExtractor>(*this,
                                                                    extractor["outputFile"].get<std::string>(),
                                                                    extractor["timeBetweenWriting"].get<double>(),
                                                                    extractor["stateToWrite"].get<unsigned short>(),
                                                                    extractor["points"].get<std::vector<std::vector<double>>>()));
        }
        else if(extractor["type"].get<std::string>() == "GMSH")
        {
            if(GMSHExtractorCount < 1)
            {
                m_pExtractor.push_back(std::make_unique<GMSHExtractor>(*this,
                                                                       extractor["outputFile"].get<std::string>(),
                                                                       extractor["timeBetweenWriting"].get<double>(),
                                                                       extractor["whatToWrite"].get<std::vector<std::string>>(),
                                                                       whatCanBeWritten,
                                                                       extractor["writeAs"].get<std::string>()));
                GMSHExtractorCount++;
            }
            else
            {
                std::cerr << "Cannot add more than one GMSH extractor!" << std::endl;
            }
        }
        else if(extractor["type"].get<std::string>() == "MinMax")
        {
            m_pExtractor.push_back(std::make_unique<MinMaxExtractor>(*this,
                                                                     extractor["outputFile"].get<std::string>(),
                                                                     extractor["timeBetweenWriting"].get<double>(),
                                                                     extractor["coordinate"].get<unsigned short>(),
                                                                     extractor["minMax"].get<std::string>()));
        }
        else if(extractor["type"].get<std::string>() == "Mass")
        {
            m_pExtractor.push_back(std::make_unique<MassExtractor>(*this,
                                                                   extractor["outputFile"].get<std::string>(),
                                                                   extractor["timeBetweenWriting"].get<double>()));
        }
        else
        {
            std::string errotText = std::string("unknown extractor type ") + extractor["Type"].get<std::string>();
            throw std::runtime_error(errotText);
        }
    }

    m_lua["rho"] = m_rho;
    m_lua["mu"] = m_mu;

    m_MPrev.resize(dim*(dim + 1), dim*(dim + 1)); m_MPrev.setZero();
    for(unsigned short k = 0 ; k < m_N.size() ; ++k)
    {
        m_MPrev += m_N[k].transpose()*m_N[k]*m_mesh.getGaussWeight(k);
    }
    m_MPrev *= m_mesh.getRefElementSize()*m_rho;

    m_FPrev.resize(dim*(dim + 1)); m_FPrev.setZero();
    for(unsigned short k = 0 ; k < m_N.size() ; ++k)
    {
        m_FPrev += m_N[k].transpose()*m_bodyForces*m_mesh.getGaussWeight(k);;
    }
    m_FPrev *= m_mesh.getRefElementSize()*m_rho;

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
}

void SolverIncompressible::applyBoundaryConditions(Eigen::VectorXd& b, const Eigen::VectorXd& qPrev)
{
    assert(m_mesh.getNodesNumber() != 0);

    const unsigned short dim = m_mesh.getDim();

    //Do not parallelize this
    for (IndexType n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        if(m_mesh.isNodeFree(n))
        {
            b(n + dim*m_mesh.getNodesNumber()) = 0;

            if(!m_mesh.isNodeBound(n))
            {
                for(unsigned short d = 0 ; d < dim - 1 ; ++d)
                {
                    b(n + d*m_mesh.getNodesNumber()) = qPrev(n + d*m_mesh.getNodesNumber());
                }

                b(n + (dim - 1)*m_mesh.getNodesNumber()) = qPrev(n + (dim - 1)*m_mesh.getNodesNumber()) - m_currentDT*m_gravity;
            }
        }

        if(m_mesh.isNodeBound(n))
        {
            std::vector<double> result;
            result = m_lua[m_mesh.getNodeType(n)](m_mesh.getNodePosition(n),
                                                  m_mesh.getNodeInitialPosition(n),
                                                  m_currentTime + m_currentDT).get<std::vector<double>>();

            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                b(n + d*m_mesh.getNodesNumber()) = result[d];
            }
        }

        if(m_strongPAtFS && m_mesh.isNodeOnFreeSurface(n))
        {
            b(n + dim*m_mesh.getNodesNumber()) = 0;
        }
    }
}

void SolverIncompressible::displaySolverParams() const noexcept
{
    std::cout << "Initial nodes number: " << m_mesh.getNodesNumber() << std::endl;
    std::cout << "Initial elements number: " << m_mesh.getElementsNumber() << std::endl;
    std::cout << "Mesh dimension: " << m_mesh.getDim() << "D" << std::endl;
    std::cout << "alpha: " << m_mesh.getAlpha() << std::endl;
    std::cout << "hchar: " << m_mesh.getHchar() << std::endl;
    std::cout << "gamma: " << m_mesh.getGamma() << std::endl;
    std::cout << "omega: " << m_mesh.getOmega() << std::endl;
    std::cout << "Eigen sparse solver: SparseLU" << std::endl;
#if defined(_OPENMP)
    std::cout << "Number of OpenMP threads: " << m_numOMPThreads << "/" << omp_get_num_procs() << std::endl;
#endif
    std::cout << "Gravity: " << m_gravity << " m/s^2" << std::endl;
    std::cout << "Density: " << m_rho << " kg/m^3" << std::endl;
    std::cout << "Viscosity: " << m_mu << " Pa s" << std::endl;
    std::cout << "Picard relative tolerance: " << m_picardRelTol  << std::endl;
    std::cout << "Maximum picard iteration number: " << m_picardMaxIter << std::endl;
    std::cout << "End simulation time: " << m_endTime << " s" << std::endl;
    std::cout << "Adapt time step: " << (m_adaptDT ? "yes" : "no") << std::endl;
    if(m_adaptDT)
    {
        std::cout << "Maximum time step: " << m_maxDT << " s" << std::endl;
        std::cout << "Time step reduction coefficient: " << m_coeffDTdecrease << std::endl;
        std::cout << "Time step increase coefficient: " << m_coeffDTincrease << std::endl;
    }
    else
        std::cout << "Time step: " << m_maxDT << " s" << std::endl;
}

void SolverIncompressible::solveProblem()
{
    std::cout   << "================================================================"
                << std::endl
                << "            EXECUTING PFEM INCOMPRESSIBLE PSPG SOLVER           "
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
        if(m_verboseOutput)
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

        if(solveCurrentTimeStep())
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
                if(m_verboseOutput)
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

bool SolverIncompressible::solveCurrentTimeStep()
{
    TimeType startTime, endTime, startTimeMeasure, endTimeMeasure;
    startTime = Clock::now();

    startTimeMeasure = Clock::now();
    const unsigned short dim = m_mesh.getDim();

    Eigen::VectorXd qPrev = getQFromNodesStates(0, m_statesNumber - 1);    //The precedent solution

    Eigen::SparseMatrix<double> A((dim+1)*m_mesh.getNodesNumber(), (dim+1)*m_mesh.getNodesNumber());    //The matrix A representing the problem: [M/dt+K -D^T; C/dt-D L].
    Eigen::VectorXd b(m_statesNumber*m_mesh.getNodesNumber());                                          //The vector b representing the problem: [M/dt*qPrev + F; H].
    Eigen::SparseMatrix<double> M(dim*m_mesh.getNodesNumber(), dim*m_mesh.getNodesNumber());            //The mass matrix.
    Eigen::SparseMatrix<double> K;                                                                      //The viscosity matrix.
    Eigen::SparseMatrix<double> D;                                                                      //The pressure matrix.
    if(m_verboseOutput)
    {
        K.resize(dim*m_mesh.getNodesNumber(), dim*m_mesh.getNodesNumber());

        D.resize(m_mesh.getNodesNumber(), dim*m_mesh.getNodesNumber());
    }
    Eigen::SparseMatrix<double> C(m_mesh.getNodesNumber(), dim*m_mesh.getNodesNumber());
    Eigen::VectorXd F(dim*m_mesh.getNodesNumber());                                                     //The volume force vector.
    Eigen::VectorXd H(m_mesh.getNodesNumber());

    std::vector<double> tauPSPG(m_mesh.getElementsNumber());                                            //tau_PSPG parameters for each element.

    m_mesh.saveNodesList();

    m_picardCurrentNumIter = 0;
    Eigen::VectorXd qIter(m_statesNumber*m_mesh.getNodesNumber()); qIter.setZero();
    Eigen::VectorXd qIterPrev(m_statesNumber*m_mesh.getNodesNumber()); qIterPrev.setZero();
    double res = std::numeric_limits<double>::max();
    endTimeMeasure = Clock::now();
    if(m_verboseOutput)
            displayDT(startTimeMeasure, endTimeMeasure, "Picard algorithm set up in ");

    while(res > m_picardRelTol)
    {
        if(m_verboseOutput)
        {
            std::cout << "Picard algorithm (mesh position) - iteration ("
                      << m_picardCurrentNumIter << ")" << std::endl;
        }

        qIterPrev = qIter;

        startTimeMeasure = Clock::now();
        computeTauPSPG(tauPSPG);
        endTimeMeasure = Clock::now();
        if(m_verboseOutput)
            displayDT(startTimeMeasure, endTimeMeasure, " * TauPSPG computed in ");

        startTimeMeasure = Clock::now();
        buildPicardSystem(A, b, M, K, D, C, F, H, tauPSPG, qPrev);
        endTimeMeasure = Clock::now();
        if(m_verboseOutput)
            displayDT(startTimeMeasure, endTimeMeasure, " * Ax = b built in ");

        startTimeMeasure = Clock::now();
        applyBoundaryConditions(b, qPrev);
        endTimeMeasure = Clock::now();
        if(m_verboseOutput)
            displayDT(startTimeMeasure, endTimeMeasure, " * BC applied in ");

        startTimeMeasure = Clock::now();
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
        endTimeMeasure = Clock::now();
        if(m_verboseOutput)
            displayDT(startTimeMeasure, endTimeMeasure, " * System solved in ");


        startTimeMeasure = Clock::now();
        setNodesStatesfromQ(qIter, 0, m_statesNumber - 1);
        Eigen::VectorXd deltaPos = qIter*m_currentDT;
        m_mesh.updateNodesPositionFromSave(std::vector<double> (deltaPos.data(), deltaPos.data() + dim*m_mesh.getNodesNumber()));
        m_picardCurrentNumIter++;
        endTimeMeasure = Clock::now();
        if(m_verboseOutput)
            displayDT(startTimeMeasure, endTimeMeasure, " * Solution updated in ");

        startTimeMeasure = Clock::now();
        if(m_picardCurrentNumIter == 1)
            res = std::numeric_limits<double>::max();
        else
        {
            double num{0};
            double den{0};

            for(IndexType n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
            {
                if(!m_mesh.isNodeFree(n))
                {
                    for(unsigned short d = 0 ; d < dim ; ++d)
                    {
                        num += (qIter(n + d*m_mesh.getNodesNumber()) - qIterPrev(n + d*m_mesh.getNodesNumber()))*(qIter(n + d*m_mesh.getNodesNumber()) - qIterPrev(n + d*m_mesh.getNodesNumber()));
                        den += qIterPrev(n + d*m_mesh.getNodesNumber())*qIterPrev(n + d*m_mesh.getNodesNumber());
                    }
                }
            }

            res = std::sqrt(num/den);
        }
        endTimeMeasure = Clock::now();
        if(m_verboseOutput)
            displayDT(startTimeMeasure, endTimeMeasure, " * Residual computed in ");

        if(m_verboseOutput)
        {
            startTimeMeasure = Clock::now();
            std::cout << " * Relative 2-norm of q: " << res << " vs "
                      << m_picardRelTol << std::endl;

            Eigen::VectorXd mom = M*(qIter.head(dim*m_mesh.getNodesNumber()) - qPrev.head(dim*m_mesh.getNodesNumber()))
                                + K*qIter.head(dim*m_mesh.getNodesNumber())
                                - Eigen::MatrixXd(D).transpose()*qIter.tail(m_mesh.getNodesNumber())
                                - F;

            Eigen::VectorXd cont = D*qIter.head(dim*m_mesh.getNodesNumber());

            for(IndexType n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
            {
                if(m_mesh.isNodeBound(n) || m_mesh.isNodeFree(n))
                {
                    for(unsigned short d = 0 ; d < dim ; ++d)
                        mom(n + d*m_mesh.getNodesNumber()) = 0;

                    cont(n) = 0;
                }
            }

            std::cout << " * Error on momentum: " << mom.norm() <<std::endl;
            std::cout << " * Error on mass: " << cont.norm() <<std::endl;
            endTimeMeasure = Clock::now();
            displayDT(startTimeMeasure, endTimeMeasure, " * Error computed in ");
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
    if(m_verboseOutput)
        displayDT(startTime, endTime, "Picard algorithm solved in ");

    startTime = Clock::now();
    m_mesh.remesh();
    endTime = Clock::now();
    if(m_verboseOutput)
        displayDT(startTime, endTime, "Remeshing done in ");

    return true;
}

void SolverIncompressible::buildPicardSystem(Eigen::SparseMatrix<double>& A,
                                             Eigen::VectorXd& b,
                                             Eigen::SparseMatrix<double>& M,
                                             Eigen::SparseMatrix<double>& K,
                                             Eigen::SparseMatrix<double>& D,
                                             Eigen::SparseMatrix<double>& C,
                                             Eigen::VectorXd& F,
                                             Eigen::VectorXd& H,
                                             const std::vector<double>& tauPSPG, const Eigen::VectorXd& qPrev)
{
    /*A = [(matrices.M)/p.dt + matrices.K, -transpose(matrices.D);...
         (matrices.C)/p.dt - matrices.D, matrices.L];

      b = [matrices.F + (matrices.M*qPrev(1:2*Nn))/p.dt;...
           matrices.H + (matrices.C*qPrev(1:2*Nn))/p.dt];
    */

    assert(tauPSPG.size() == m_mesh.getElementsNumber());

    const unsigned short dim = m_mesh.getDim();
    const unsigned short noPerEl = dim + 1;
    const unsigned int tripletPerElm = (dim*noPerEl*noPerEl + dim*noPerEl*dim*noPerEl + 3*noPerEl*dim*noPerEl + noPerEl*noPerEl);
    const IndexType nElm = m_mesh.getElementsNumber();
    const IndexType nNodes = m_mesh.getNodesNumber();

    std::vector<Eigen::Triplet<double>> indexA(tripletPerElm*nElm);

    std::vector<Eigen::Triplet<double>> indexM(dim*noPerEl*noPerEl*nElm);

    std::vector<Eigen::Triplet<double>> indexK;
    std::vector<Eigen::Triplet<double>> indexD;
    if(m_verboseOutput)
    {
        indexK.resize(dim*noPerEl*dim*noPerEl*nElm);

        indexD.resize(dim*noPerEl*noPerEl*nElm);
    }

    std::vector<Eigen::Triplet<double>> indexC(dim*noPerEl*noPerEl*nElm);

    std::vector<std::pair<IndexType, double>> indexF(dim*(dim + 1)*nElm, std::make_pair(0, 0.0));
    std::vector<std::pair<IndexType, double>> indexH((dim + 1)*nElm, std::make_pair(0, 0.0));

    F.setZero();
    H.setZero();

    Eigen::MatrixXd MPrev = m_MPrev/m_currentDT;

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(IndexType elm = 0 ; elm < nElm ; ++elm)
    {
        Eigen::MatrixXd Be = getB(elm);
        Eigen::MatrixXd Bep(dim, dim + 1);
        if(dim == 2)
            Bep << Be.block(0, 0, 1, 3), Be.block(1, 3, 1, 3);
        else
            Bep << Be.block(0, 0, 1, 4), Be.block(1, 4, 1, 4), Be.block(2, 8, 1, 4);

        //Me = S rho Nv^T Nv dV
        Eigen::MatrixXd Me = MPrev*m_mesh.getElementDetJ(elm);

        //Ke = S Bv^T ddev Bv dV
        Eigen::MatrixXd Ke = m_mesh.getRefElementSize()*Be.transpose()*m_ddev*Be*m_mesh.getElementDetJ(elm); //same matrices for all gauss point ^^

        //De = S Np^T m Bv dV
        Eigen::MatrixXd De((dim+1),dim*(dim+1)); De.setZero();

        for (unsigned short k = 0 ; k < m_N.size() ; ++k)
            De += (m_N[k].topLeftCorner(1, dim+1)).transpose()*m_m.transpose()*Be*m_mesh.getGaussWeight(k);


        De *= m_mesh.getRefElementSize()*m_mesh.getElementDetJ(elm);

        //Ce = tauPSPG S Bp^T Nv dV
        Eigen::MatrixXd Ce((dim+1),dim*(dim+1)); Ce.setZero();

        for (unsigned short k = 0 ; k < m_N.size() ; ++k)
            Ce += Bep.transpose()*m_N[k]*m_mesh.getGaussWeight(k);

        Ce *= (m_mesh.getRefElementSize()*tauPSPG[elm]/m_currentDT)*m_mesh.getElementDetJ(elm);

        //Le = tauPSPG S Bp^T Bp dV
        Eigen::MatrixXd Le = m_mesh.getRefElementSize()*(tauPSPG[elm]/m_rho)*
                             Bep.transpose()*Bep*m_mesh.getElementDetJ(elm);

        //Fe = S Nv^T bodyforce dV
        Eigen::MatrixXd Fe = m_FPrev*m_mesh.getElementDetJ(elm);

        //He = tauPSPG S Bp^T bodyforce dV
        Eigen::MatrixXd He = m_mesh.getRefElementSize()*tauPSPG[elm]*Bep.transpose()*
                             m_bodyForces*m_mesh.getElementDetJ(elm);

        std::size_t countA = 0;
        std::size_t countM = 0;
        std::size_t countCD = 0;
        std::size_t countK = 0;
        std::size_t countH = 0;
        std::size_t countF = 0;

        for(unsigned short i = 0 ; i < (dim+1) ; ++i)
        {
            for(unsigned short j = 0 ; j < (dim+1) ; ++j)
            {
                for(unsigned short d = 0 ; d < dim ; ++d)
                {
                    /********************************************************************
                                                 Build M/dt
                    ********************************************************************/
                    indexM[dim*noPerEl*noPerEl*elm + countM] =
                        Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + d*nNodes,
                                               m_mesh.getElement(elm)[j] + d*nNodes,
                                               Me(i + d*noPerEl, j + d*noPerEl));

                    countM++;

                    if(!(m_mesh.isNodeBound(m_mesh.getElement(elm)[i]) || m_mesh.isNodeFree(m_mesh.getElement(elm)[i])))
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + d*nNodes,
                                                   m_mesh.getElement(elm)[j] + d*nNodes,
                                                   Me(i + d*noPerEl, j + d*noPerEl));
                    }

                    countA++;

                    /********************************************************************
                                                  Build K
                    ********************************************************************/
                    for(unsigned short d2 = 0 ; d2 < dim ; ++d2)
                    {
                        if(m_verboseOutput)
                        {
                            indexK[dim*noPerEl*dim*noPerEl*elm + countK] =
                                Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + d*nNodes,
                                                       m_mesh.getElement(elm)[j] + d2*nNodes,
                                                       Ke(i + d*noPerEl, j + d2*noPerEl));

                            countK++;
                        }

                        if(!(m_mesh.isNodeBound(m_mesh.getElement(elm)[i]) || m_mesh.isNodeFree(m_mesh.getElement(elm)[i])))
                        {
                            indexA[tripletPerElm*elm + countA] =
                                Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + d*nNodes,
                                                       m_mesh.getElement(elm)[j] + d2*nNodes,
                                                       Ke(i + d*noPerEl, j + d2*noPerEl));
                        }
                        countA++;
                    }

                    /********************************************************************
                                                  Build D
                    ********************************************************************/
                    if(m_verboseOutput)
                    {
                        indexD[noPerEl*dim*noPerEl*elm + countCD] =
                            Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                   m_mesh.getElement(elm)[j] + d*nNodes,
                                                   De(i, j + d*noPerEl));

                        //Updated with C
                    }

                    if(!(m_mesh.isNodeFree(m_mesh.getElement(elm)[i]) || (m_strongPAtFS && m_mesh.isNodeOnFreeSurface(m_mesh.getElement(elm)[i]))))
                    {
                        //Part D of A
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + dim*nNodes,
                                                   m_mesh.getElement(elm)[j] + d*nNodes,
                                                   De(i, j + d*noPerEl));
                    }

                    countA++;

                    //Part -D^T of A
                    if(!(m_mesh.isNodeBound(m_mesh.getElement(elm)[j]) || m_mesh.isNodeFree(m_mesh.getElement(elm)[j])))
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(m_mesh.getElement(elm)[j] + d*nNodes,
                                                   m_mesh.getElement(elm)[i] + dim*nNodes,
                                                   -De(i, j + d*noPerEl));
                    }

                    countA++;

                    /********************************************************************
                                                Build C/dt
                    ********************************************************************/
                    indexC[noPerEl*dim*noPerEl*elm + countCD] =
                        Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                               m_mesh.getElement(elm)[j] + d*nNodes,
                                               Ce(i, j + d*noPerEl));

                    countCD++;

                    if(!(m_mesh.isNodeFree(m_mesh.getElement(elm)[i]) || (m_strongPAtFS && m_mesh.isNodeOnFreeSurface(m_mesh.getElement(elm)[i]))))
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + dim*nNodes,
                                                   m_mesh.getElement(elm)[j] + d*nNodes,
                                                   Ce(i, j + d*noPerEl));
                    }

                    countA++;
                }

                /********************************************************************
                                            Build L
                ********************************************************************/
                if(!(m_mesh.isNodeFree(m_mesh.getElement(elm)[i]) || (m_strongPAtFS && m_mesh.isNodeOnFreeSurface(m_mesh.getElement(elm)[i]))))
                {
                    indexA[tripletPerElm*elm + countA] =
                        Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + dim*nNodes,
                                               m_mesh.getElement(elm)[j] + dim*nNodes,
                                               Le(i,j));
                }

                countA++;
            }

            /************************************************************************
                                              Build f
            ************************************************************************/
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                indexF[dim*(dim + 1)*elm + countF] =
                    std::make_pair(m_mesh.getElement(elm)[i] + d*nNodes, Fe(i + d*noPerEl));

                countF++;
            }

            /************************************************************************
                                            Build h
            ************************************************************************/
            indexH[(dim + 1)*elm + countH] =
                    std::make_pair(m_mesh.getElement(elm)[i], He(i));

            countH++;
        }
    }
    Eigen::setNbThreads(m_numOMPThreads);

    //Best would be to know the number of nodes in which case :/
    //This can still be fasten using openmp but will never be as good as using []
    //with preallocated memory
    for(IndexType n = 0 ; n < nNodes ; ++n)
    {
        if(m_mesh.isNodeFree(n) || (m_strongPAtFS && m_mesh.isNodeOnFreeSurface(n)))
        {
            indexA.push_back(Eigen::Triplet<double>(n + dim*nNodes,
                                                    n + dim*nNodes,
                                                    1));
        }

        if(m_mesh.isNodeBound(n) || m_mesh.isNodeFree(n))
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                indexA.push_back(Eigen::Triplet<double>(n + d*nNodes,
                                                        n + d*nNodes,
                                                        1));
            }
        }
    }

    for(auto& doublet : indexF)
    {
        F[doublet.first] += doublet.second;
    }

    for(auto& doublet : indexH)
    {
        H[doublet.first] += doublet.second;
    }

    M.setFromTriplets(indexM.begin(), indexM.end());

    if(m_verboseOutput)
    {
        K.setFromTriplets(indexK.begin(), indexK.end());
        D.setFromTriplets(indexD.begin(), indexD.end());
    }

    C.setFromTriplets(indexC.begin(), indexC.end());

    /********************************************************************************
                                        Compute A and b
    ********************************************************************************/
    A.setFromTriplets(indexA.begin(), indexA.end());

    b << F + M*qPrev.head(dim*nNodes), H + C*qPrev.head(dim*nNodes);
}

//Put this is build PicardSystem ?
void SolverIncompressible::computeTauPSPG(std::vector<double>& tauPSPG)
{
    #pragma omp parallel for default(shared)
    for(IndexType elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        const double h = std::sqrt(m_mesh.getRefElementSize()*m_mesh.getElementDetJ(elm)/M_PI);

        double U = 0;
        for (unsigned short n = 0 ; n < m_mesh.getDim() + 1 ; ++n)
        {
            double nodeU = 0;
            for (unsigned short d = 0 ; d < m_mesh.getDim() ; ++d)
            {
                nodeU += m_mesh.getNodeState(m_mesh.getElement(elm)[n], d)*m_mesh.getNodeState(m_mesh.getElement(elm)[n], d);
            }
            U += std::sqrt(nodeU);
        }
        U /= (m_mesh.getDim() + 1);

        tauPSPG[elm] = 1/std::sqrt((2/m_currentDT)*(2/m_currentDT)
                                 + (2*U/h)*(2*U/h)
                                 + 9*(4*m_mu/(h*h*m_rho))*(4*m_mu/(h*h*m_rho)));
    }
}
