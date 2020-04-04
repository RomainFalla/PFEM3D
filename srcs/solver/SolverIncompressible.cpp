#include <iomanip>
#include <iostream>

#include <gmsh.h>

#include "SolverIncompressible.hpp"


SolverIncompressible::SolverIncompressible(const nlohmann::json& j, const std::string& mshName, const std::string& resultsName) :
Solver(j, mshName, resultsName)
{
    unsigned short dim = m_mesh.getMeshDim();

    m_statesNumber = dim + 1;

    m_mesh.setStatesNumber(m_statesNumber);

    if(m_initialCondition.size() != m_statesNumber)
        throw std::runtime_error("Invalid number of initial condition!");

    m_rho               = j["Solver"]["Fluid"]["rho"].get<double>();
    m_mu                = j["Solver"]["Fluid"]["mu"].get<double>();

    m_picardRelTol           = j["Solver"]["Picard"]["relTol"].get<double>();
    m_picardMaxIter          = j["Solver"]["Picard"]["maxIter"].get<unsigned int>();

    m_coeffDTincrease    = j["Solver"]["Time"]["coeffDTincrease"].get<double>();
    m_coeffDTdecrease    = j["Solver"]["Time"]["coeffDTdecrease"].get<double>();

    if(dim == 2)
        m_whatCanBeWriten = {"u", "v", "p", "ke", "velocity"};
    else if(dim == 3)
        m_whatCanBeWriten = {"u", "v", "w", "p", "ke", "velocity"};

    std::vector<std::string> whatToWrite = j["Solver"]["whatToWrite"];
    m_whatToWrite.resize(dim + 3, false);

    for(unsigned short i = 0 ; i < whatToWrite.size() ; ++i)
    {
        bool found = false;
        for(unsigned short j = 0 ; j < m_whatCanBeWriten.size() ; ++j)
        {
            if(whatToWrite[i] == m_whatCanBeWriten[j])
            {
                m_whatToWrite[j] = true;
                found = true;
                break;
            }
        }
        if(!found)
            throw std::runtime_error("Unknown quantity to write!");
    }

    gmsh::initialize();

#ifndef NDEBUG
    gmsh::option::setNumber("General.Terminal", 1);
#else
    gmsh::option::setNumber("General.Terminal", 0);
#endif // DEBUG

    gmsh::option::setNumber("General.NumThreads", m_numOMPThreads);

    m_sumNTN.resize(dim*(dim + 1), dim*(dim + 1)); m_sumNTN.setZero();
    for(unsigned short k = 0 ; k < m_N.size() ; ++k)
    {
        m_sumNTN += m_N[k].transpose()*m_N[k]*m_mesh.getGaussWeight(k);
    }

    if(dim == 2)
    {
        m_ddev.resize(3, 3);
        m_ddev << 2, 0, 0,
                  0, 2, 0,
                  0, 0, 1;
    }
    else if(dim == 3)
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
}

SolverIncompressible::~SolverIncompressible()
{
    gmsh::finalize();
}

void SolverIncompressible::applyBoundaryConditions()
{
    assert(m_mesh.getNodesNumber() != 0);

    const unsigned short dim = m_mesh.getMeshDim();

    //Do not parallelize this
    for (std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        if(m_mesh.isNodeFree(n))
        {
            m_A.row(n + dim*m_mesh.getNodesNumber()) *= 0;
            m_b(n + dim*m_mesh.getNodesNumber()) = m_qprev(n + dim*m_mesh.getNodesNumber());
            m_A.coeffRef(n + dim*m_mesh.getNodesNumber(), n + dim*m_mesh.getNodesNumber()) = 1;

            if(!m_mesh.isNodeBound(n))
            {
                for(unsigned short d = 0 ; d < dim - 1 ; ++d)
                {
                    m_A.row(n + d*m_mesh.getNodesNumber()) *= 0;
                    m_b(n + d*m_mesh.getNodesNumber()) = m_qprev(n + d*m_mesh.getNodesNumber());
                    m_A.coeffRef(n + d*m_mesh.getNodesNumber(), n + d*m_mesh.getNodesNumber()) = 1;
                }

                m_A.row(n + (dim - 1)*m_mesh.getNodesNumber()) *= 0;
                m_b(n + (dim - 1)*m_mesh.getNodesNumber()) = m_qprev(n + (dim - 1)*m_mesh.getNodesNumber()) - m_currentDT*m_gravity;
                m_A.coeffRef(n + (dim - 1)*m_mesh.getNodesNumber(), n + (dim - 1)*m_mesh.getNodesNumber()) = 1;
            }
        }

        if(m_mesh.isNodeBound(n))
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                m_A.row(n + d*m_mesh.getNodesNumber()) *= 0;
                m_b(n + d*m_mesh.getNodesNumber()) = m_qprev(n + d*m_mesh.getNodesNumber());
                m_A.coeffRef(n + d*m_mesh.getNodesNumber(), n + d*m_mesh.getNodesNumber()) = 1;
            }
        }
    }

    m_A.prune(0, 0);
}

void SolverIncompressible::displaySolverParams() const
{
    std::cout << "Initial nodes number: " << m_mesh.getNodesNumber() << std::endl;
    std::cout << "Initial elements number: " << m_mesh.getElementsNumber() << std::endl;
    std::cout << "Mesh dimension: " << m_mesh.getMeshDim() << "D" << std::endl;
    std::cout << "alpha: " << m_mesh.getAlpha() << std::endl;
    std::cout << "hchar: " << m_mesh.getHchar() << std::endl;
    std::cout << "gamma: " << m_mesh.getGamma() << std::endl;
    std::cout << "omega: " << m_mesh.getOmega() << std::endl;
    std::cout << "Eigen sparse solver: SparseLU" << std::endl;
#if defined(_OPENMP)
    std::cout << "Number of OpenMP threads: " << m_numOMPThreads << std::endl;
#endif
    std::cout << "Gravity: " << m_gravity << " m/s^2" << std::endl;
    std::cout << "Density: " << m_rho << " kg/m^3" << std::endl;
    std::cout << "Viscosity: " << m_mu << " Pa s" << std::endl;
    std::cout << "Picard relative tolerance: " << m_picardRelTol  << std::endl;
    std::cout << "Maximum picard iteration number: " << m_picardMaxIter << std::endl;
    std::cout << "End simulation time: " << m_endTime << " s" << std::endl;
    std::cout << "Write solution every: " << m_timeBetweenWriting << " s" << std::endl;
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

void SolverIncompressible::setInitialCondition()
{
    assert(m_mesh.getNodesNumber() != 0);

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        for(unsigned short s = 0 ; s < m_mesh.getMeshDim() + 1 ; ++s)
        {
            if(!m_mesh.isNodeBound(n) || m_mesh.isNodeFluidInput(n))
            {
                m_mesh.setNodeState(n, s, m_initialCondition[s]);
            }
            else
            {
                m_mesh.setNodeState(n, s, 0);
            }
        }
    }

    m_qprev = getQFromNodesStates(0, m_statesNumber - 1);
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

    setInitialCondition();

    writeData();

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
            if(m_currentTime >= m_nextWriteTrigger)
            {
                writeData();
                m_nextWriteTrigger +=m_timeBetweenWriting;
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
    assert(m_qprev.size() == m_statesNumber*m_mesh.getNodesNumber());

    const unsigned short dim = m_mesh.getMeshDim();

    m_mesh.saveNodesList();

    m_picardCurrentNumIter = 0;
    Eigen::VectorXd qIter(m_statesNumber*m_mesh.getNodesNumber()); qIter.setZero();
    Eigen::VectorXd qIterPrev(m_statesNumber*m_mesh.getNodesNumber()); qIterPrev.setZero();
    double res = std::numeric_limits<double>::max();

    while(res > m_picardRelTol)
    {
        if(m_verboseOutput)
        {
            std::cout << "Picard algorithm (mesh position) - iteration ("
                      << m_picardCurrentNumIter << ")" << std::endl;
        }

        qIterPrev = qIter;

        computeTauPSPG();
        buildPicardSystem();
        applyBoundaryConditions();

        m_solverLU.compute(m_A);
        if(m_solverLU.info() == Eigen::Success)
        {
            qIter = m_solverLU.solve(m_b);
        }
        else
        {
            m_mesh.restoreNodesList();
            return false;
        }

        setNodesStatesfromQ(qIter, 0, m_statesNumber - 1);
        Eigen::VectorXd deltaPos = qIter*m_currentDT;
        m_mesh.updateNodesPositionFromSave(std::vector<double> (deltaPos.data(), deltaPos.data() + dim*m_mesh.getNodesNumber()));
        m_picardCurrentNumIter++;

        if(m_picardCurrentNumIter == 1)
            res = std::numeric_limits<double>::max();
        else
        {
            double num{0};
            double den{0};

            for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
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

        if(m_verboseOutput)
        {
            std::cout << " * Relative 2-norm of q: " << res << " vs "
                      << m_picardRelTol << std::endl;

            Eigen::VectorXd mom = m_M*(qIter.head(dim*m_mesh.getNodesNumber()) - m_qprev.head(dim*m_mesh.getNodesNumber()))
                                + m_K*qIter.head(dim*m_mesh.getNodesNumber())
                                - Eigen::MatrixXd(m_D).transpose()*qIter.tail(m_mesh.getNodesNumber())
                                - m_F;

            Eigen::VectorXd cont = m_D*qIter.head(dim*m_mesh.getNodesNumber());

            for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
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
        }

        if(m_picardCurrentNumIter >= m_picardMaxIter || std::isnan(res))
        {
            m_mesh.restoreNodesList();
            return false;
        }
    }

    m_currentTime += m_currentDT;
    m_currentStep++;

    //Remeshing step
    m_mesh.remesh();

    //We have to compute qPrev here due to new nodes !
    m_qprev = getQFromNodesStates(0, m_statesNumber - 1);

    return true;
}

void SolverIncompressible::buildPicardSystem()
{
    assert(m_tauPSPG.size() == m_mesh.getElementsNumber());

    const unsigned short dim = m_mesh.getMeshDim();
    const unsigned short noPerEl = dim + 1;

    /*A = [(matrices.M)/p.dt + matrices.K, -transpose(matrices.D);...
         (matrices.C)/p.dt - matrices.D, matrices.L];

      b = [matrices.F + (matrices.M*qPrev(1:2*Nn))/p.dt;...
           matrices.H + (matrices.C*qPrev(1:2*Nn))/p.dt];
    */

    m_A.resize((dim+1)*m_mesh.getNodesNumber(), (dim+1)*m_mesh.getNodesNumber());
    m_A.data().squeeze();
    std::vector<Eigen::Triplet<double>> indexA;
    indexA.reserve((dim*noPerEl*noPerEl + dim*noPerEl*dim*noPerEl + 3*noPerEl*dim*noPerEl + noPerEl*noPerEl)*m_mesh.getElementsNumber());

    m_b.resize(m_statesNumber*m_mesh.getNodesNumber());
    m_b.setZero();

    m_M.resize(dim*m_mesh.getNodesNumber(), dim*m_mesh.getNodesNumber());
    m_M.data().squeeze();
    std::vector<Eigen::Triplet<double>> indexM;
    indexM.reserve(dim*noPerEl*noPerEl*m_mesh.getElementsNumber());

    std::vector<Eigen::Triplet<double>> indexK;
    std::vector<Eigen::Triplet<double>> indexD;
    if(m_verboseOutput)
    {
        m_K.resize(dim*m_mesh.getNodesNumber(), dim*m_mesh.getNodesNumber());
        m_K.data().squeeze();
        indexK.reserve(dim*noPerEl*dim*noPerEl*m_mesh.getElementsNumber());

        m_D.resize(m_mesh.getNodesNumber(), dim*m_mesh.getNodesNumber());
        m_D.data().squeeze();
        indexD.reserve(dim*noPerEl*noPerEl*m_mesh.getElementsNumber());
    }

    Eigen::SparseMatrix<double> C(m_mesh.getNodesNumber(), dim*m_mesh.getNodesNumber());
    std::vector<Eigen::Triplet<double>> indexC;
    indexC.reserve(dim*noPerEl*noPerEl*m_mesh.getElementsNumber());

    m_F.resize(dim*m_mesh.getNodesNumber());
    m_F.setZero();

    Eigen::VectorXd H(m_mesh.getNodesNumber());
    H.setZero();

    Eigen::MatrixXd MPrev = m_sumNTN*m_mesh.getRefElementSize()*(m_rho/m_currentDT);

    std::vector<Eigen::MatrixXd> Me(m_mesh.getElementsNumber());
    std::vector<Eigen::MatrixXd> Ke(m_mesh.getElementsNumber());
    std::vector<Eigen::MatrixXd> De(m_mesh.getElementsNumber());
    std::vector<Eigen::MatrixXd> Ce(m_mesh.getElementsNumber());
    std::vector<Eigen::MatrixXd> Le(m_mesh.getElementsNumber());
    std::vector<Eigen::MatrixXd> Fe(m_mesh.getElementsNumber());
    std::vector<Eigen::MatrixXd> He(m_mesh.getElementsNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        Eigen::MatrixXd Be = getB(elm);
        Eigen::MatrixXd Bep(dim, dim + 1);
        if(dim == 2)
            Bep << Be.block(0, 0, 1, 3), Be.block(1, 3, 1, 3);
        else if(dim == 3)
            Bep << Be.block(0, 0, 1, 4), Be.block(1, 4, 1, 4), Be.block(2, 8, 1, 4);

        //Me = S rho Nv^T Nv dV
        Me[elm] = MPrev*m_mesh.getElementDetJ(elm);

        //Ke = S Bv^T ddev Bv dV
        Ke[elm] = m_mesh.getRefElementSize()*Be.transpose()*m_ddev*Be*m_mesh.getElementDetJ(elm); //same matrices for all gauss point ^^

        //De = S Np^T m Bv dV
        De[elm].resize((dim+1),dim*(dim+1)); De[elm].setZero();

        for (unsigned short k = 0 ; k < m_N.size() ; ++k)
            De[elm] += (m_N[k].topLeftCorner(1, dim+1)).transpose()*m_m.transpose()*Be*m_mesh.getGaussWeight(k);


        De[elm] *= m_mesh.getRefElementSize()*m_mesh.getElementDetJ(elm);

        //Ce = tauPSPG S Bp^T Nv dV
        Ce[elm].resize((dim+1),dim*(dim+1)); Ce[elm].setZero();

        for (unsigned short k = 0 ; k < m_N.size() ; ++k)
            Ce[elm] += Bep.transpose()*m_N[k]*m_mesh.getGaussWeight(k);

        Ce[elm] *= (m_mesh.getRefElementSize()*m_tauPSPG[elm]/m_currentDT)*m_mesh.getElementDetJ(elm);

        //Le = tauPSPG S Bp^T Bp dV
        Le[elm] = m_mesh.getRefElementSize()*(m_tauPSPG[elm]/m_rho)
                *Bep.transpose()*Bep
                *m_mesh.getElementDetJ(elm);

        //Fe = S Nv^T bodyforce dV
        Fe[elm].resize(dim*(dim+1),1); Fe[elm].setZero();

        for (unsigned short k = 0 ; k < m_N.size() ; ++k)
            Fe[elm] += m_N[k].transpose()*m_bodyForces*m_mesh.getGaussWeight(k);

        Fe[elm] *= m_mesh.getRefElementSize()*m_rho*m_mesh.getElementDetJ(elm);

        //He = tauPSPG S Bp^T bodyforce dV
        He[elm] = m_mesh.getRefElementSize()*m_tauPSPG[elm]*Bep.transpose()
                *m_bodyForces*m_mesh.getElementDetJ(elm);
    }

    //Big push_back ^^
    for(std::size_t elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        for(unsigned short i = 0 ; i < (dim+1) ; ++i)
        {
            for(unsigned short j = 0 ; j < (dim+1) ; ++j)
            {
                for(unsigned short d = 0 ; d < dim ; ++d)
                {
                    /********************************************************************
                                                 Build M/dt
                    ********************************************************************/
                    indexM.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + d*m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j] + d*m_mesh.getNodesNumber(),
                                                            Me[elm](i + d*noPerEl, j + d*noPerEl)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + d*m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j] + d*m_mesh.getNodesNumber(),
                                                            Me[elm](i + d*noPerEl, j + d*noPerEl)));

                    /********************************************************************
                                                  Build K
                    ********************************************************************/
                    for(unsigned short d2 = 0 ; d2 < dim ; ++d2)
                    {
                        if(m_verboseOutput)
                        {
                            indexK.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + d*m_mesh.getNodesNumber(),
                                                                    m_mesh.getElement(elm)[j] + d2*m_mesh.getNodesNumber(),
                                                                    Ke[elm](i + d*noPerEl, j + d2*noPerEl)));
                        }

                        indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + d*m_mesh.getNodesNumber(),
                                                                m_mesh.getElement(elm)[j] + d2*m_mesh.getNodesNumber(),
                                                                Ke[elm](i + d*noPerEl, j + d2*noPerEl)));
                    }

                    /********************************************************************
                                                  Build D
                    ********************************************************************/
                    if(m_verboseOutput)
                    {
                        indexD.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                                m_mesh.getElement(elm)[j] + d*m_mesh.getNodesNumber(),
                                                                De[elm](i, j + d*noPerEl)));
                    }

                    //Part D of A
                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + dim*m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j] + d*m_mesh.getNodesNumber(),
                                                            De[elm](i, j + d*noPerEl)));

                    //Part -D^T of A
                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[j] + d*m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[i] + dim*m_mesh.getNodesNumber(),
                                                            -De[elm](i, j + d*noPerEl)));

                    /********************************************************************
                                                Build C/dt
                    ********************************************************************/
                    indexC.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j] + d*m_mesh.getNodesNumber(),
                                                            Ce[elm](i, j + d*noPerEl)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + dim*m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j] + d*m_mesh.getNodesNumber(),
                                                            Ce[elm](i, j + d*noPerEl)));
                }

                /********************************************************************
                                            Build L
                ********************************************************************/
                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + dim*m_mesh.getNodesNumber(),
                                                        m_mesh.getElement(elm)[j] + dim*m_mesh.getNodesNumber(),
                                                        Le[elm](i,j)));
            }

            /************************************************************************
                                              Build f
            ************************************************************************/
            for(unsigned short d = 0 ; d < dim ; ++d)
                m_F(m_mesh.getElement(elm)[i] + d*m_mesh.getNodesNumber()) += Fe[elm](i + d*noPerEl);

            /************************************************************************
                                            Build h
            ************************************************************************/
            H(m_mesh.getElement(elm)[i]) += He[elm](i);
        }
    }

    m_M.setFromTriplets(indexM.begin(), indexM.end());

    if(m_verboseOutput)
    {
        m_K.setFromTriplets(indexK.begin(), indexK.end());
        m_D.setFromTriplets(indexD.begin(), indexD.end());
    }

    C.setFromTriplets(indexC.begin(), indexC.end());

    /********************************************************************************
                                        Compute A and b
    ********************************************************************************/
    m_A.setFromTriplets(indexA.begin(), indexA.end());

    m_b << m_F + m_M*m_qprev.head(dim*m_mesh.getNodesNumber()), H + C*m_qprev.head(dim*m_mesh.getNodesNumber());
}

void SolverIncompressible::computeTauPSPG()
{
    m_tauPSPG.resize(m_mesh.getElementsNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        const double h = std::sqrt(m_mesh.getRefElementSize()*m_mesh.getElementDetJ(elm)/M_PI);

        double U = 0;
        for (unsigned short n = 0 ; n < m_mesh.getMeshDim() + 1 ; ++n)
        {
            double nodeU = 0;
            for (unsigned short d = 0 ; d < m_mesh.getMeshDim() ; ++d)
            {
                nodeU += m_mesh.getNodeState(m_mesh.getElement(elm)[n], d)*m_mesh.getNodeState(m_mesh.getElement(elm)[n], d);
            }
            U += std::sqrt(nodeU);
        }
        U /= (m_mesh.getMeshDim() + 1);

        m_tauPSPG[elm] = 1/std::sqrt((2/m_currentDT)*(2/m_currentDT)
                                 + (2*U/h)*(2*U/h)
                                 + 9*(4*m_mu/(h*h*m_rho))*(4*m_mu/(h*h*m_rho)));
    }
}
