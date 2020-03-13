#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#if defined(_OPENMP)
    #include <cstdlib>
    #include <omp.h>
#endif

#include <gmsh.h>
#include <nlohmann/json.hpp>

#include "Solver.hpp"
#include "../quadrature/gausslegendre.hpp"


Solver::Solver(const nlohmann::json& j, std::string mshName, std::string resultsName) :
m_mesh(j)
{
    m_verboseOutput             = j["verboseOutput"].get<bool>();

    m_p.hchar                   = j["RemeshingParams"]["hchar"].get<double>();

    m_p.gravity                 = j["SolverIncompressibleParams"]["gravity"].get<double>();

    m_p.fluid.rho               = j["SolverIncompressibleParams"]["FluidParams"]["rho"].get<double>();
    m_p.fluid.mu                = j["SolverIncompressibleParams"]["FluidParams"]["mu"].get<double>();

    m_p.picard.relTol           = j["SolverIncompressibleParams"]["PicardParams"]["relTol"].get<double>();
    m_p.picard.maxIter          = j["SolverIncompressibleParams"]["PicardParams"]["maxIter"].get<unsigned int>();

    m_p.time.adaptDT            = j["SolverIncompressibleParams"]["TimeParams"]["adaptDT"].get<bool>();
    m_p.time.coeffDTincrease    = j["SolverIncompressibleParams"]["TimeParams"]["coeffDTincrease"].get<double>();
    m_p.time.coeffDTdecrease    = j["SolverIncompressibleParams"]["TimeParams"]["coeffDTdecrease"].get<double>();
    m_p.time.maxDT              = j["SolverIncompressibleParams"]["TimeParams"]["maxDT"].get<double>();
    m_p.time.simuTime           = j["SolverIncompressibleParams"]["TimeParams"]["simuTime"].get<double>();
    m_p.time.simuDTToWrite      = j["SolverIncompressibleParams"]["TimeParams"]["simuDTToWrite"].get<double>();
    m_p.time.currentDT          = m_p.time.maxDT;
    m_p.time.nextWriteTrigger   = m_p.time.simuDTToWrite;
    m_p.time.currentTime        = 0.0;
    m_p.time.currentStep        = 0;

    m_p.initialCondition        = j["SolverIncompressibleParams"]["initialCondition"].get<std::vector<double>>();

    // set the desired number of OpenMP threads
#if defined(_OPENMP)
    const char* pNumThreads = std::getenv("OMP_NUM_THREADS");

    if(pNumThreads == nullptr)
        m_p.nOMPThreads = 1;
    else
        m_p.nOMPThreads = std::atoi(pNumThreads);

    omp_set_num_threads(m_p.nOMPThreads);
    Eigen::setNbThreads(m_p.nOMPThreads);
#else
     m_p.nOMPThreads = 1;
#endif

    m_mesh.loadFromFile(mshName);

    std::vector<std::string> whatToWrite = j["SolverIncompressibleParams"]["whatToWrite"];
    for(auto what : whatToWrite)
    {
        if(what == "u")
            m_p.whatToWrite[0] = true;
        else if(what == "v")
            m_p.whatToWrite[1] = true;
        else if(what == "p")
            m_p.whatToWrite[2] = true;
        else if(what == "ke")
            m_p.whatToWrite[3] = true;
        else if(what == "velocity")
            m_p.whatToWrite[4] = true;
        else
            throw std::runtime_error("Unknown quantity to write!");
    }

    m_p.writeAs = j["SolverIncompressibleParams"]["writeAs"];
    if(!(m_p.writeAs == "Nodes" || m_p.writeAs == "Elements" || m_p.writeAs == "NodesElements"))
        throw std::runtime_error("Unexpected date type to write!");

    gmsh::initialize();

#ifndef NDEBUG
    gmsh::option::setNumber("General.Terminal", 1);
#else
    gmsh::option::setNumber("General.Terminal", 0);
#endif // DEBUG

    gmsh::option::setNumber("General.NumThreads", m_p.nOMPThreads);

    m_p.resultsName = resultsName;

    m_N = getN();

    m_sumNTN.resize(6,6); m_sumNTN.setZero();
    for(unsigned short k = 0 ; k < m_N.size() ; ++k)
        m_sumNTN += m_N[k].transpose()*m_N[k]*GP2Dweight<double>[k];

    m_m.resize(3);
    m_m << 1, 1, 0;

    m_ddev.resize(3, 3);
    m_ddev << 2, 0, 0,
              0, 2, 0,
              0, 0, 1;

    m_ddev *= m_p.fluid.mu;
}

Solver::~Solver()
{
    gmsh::finalize();
}

void Solver::applyBoundaryConditions()
{
	  assert(!m_mesh.nodesList.empty());

    m_A.prune([this](int i, int j, float)
    {
        if(i < m_mesh.getNodesNumber())
        {
            return !(this->m_mesh.isNodeBound(i) || this->m_mesh.isNodeFree(i));
        }
        else if (i >= m_mesh.getNodesNumber() && i < 2*m_mesh.getNodesNumber())
        {
            return !(this->m_mesh.isNodeBound(i - m_mesh.getNodesNumber()) ||
                    this->m_mesh.isNodeFree(i - m_mesh.getNodesNumber()));
        }
        else
        {
            return !(this->m_mesh.isNodeFree(i - 2*m_mesh.getNodesNumber()));
        }
    });

    //Do not parallelize this
    for (std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        if(m_mesh.isNodeFree(n))
        {
            m_b(n + 2*m_mesh.getNodesNumber()) = m_qprev(n + 2*m_mesh.getNodesNumber());
            m_A.coeffRef(n + 2*m_mesh.getNodesNumber(), n + 2*m_mesh.getNodesNumber()) = 1;

            if(!m_mesh.isNodeBound(n))
            {
                m_b(n) = m_qprev(n);
                m_A.coeffRef(n, n) = 1;

                m_b(n + m_mesh.getNodesNumber()) = m_qprev(n + m_mesh.getNodesNumber()) - m_p.time.currentDT*m_p.gravity*m_p.fluid.rho*m_p.hchar*m_p.hchar*0.5;;
                m_A.coeffRef(n + m_mesh.getNodesNumber(), n + m_mesh.getNodesNumber()) = 1;
            }
        }

        if(m_mesh.isNodeBound(n))
        {
            m_b(n) = m_qprev(n);
            m_A.coeffRef(n, n) = 1;

            m_b(n + m_mesh.getNodesNumber()) = m_qprev(n + m_mesh.getNodesNumber());
            m_A.coeffRef(n + m_mesh.getNodesNumber(), n + m_mesh.getNodesNumber()) = 1;
        }
    }
}

void Solver::displaySolverParams() const
{
    std::cout << "Eigen sparse solver: SparseLU" << std::endl;
#if defined(_OPENMP)
    std::cout << "Number of OpenMP threads: " << m_p.nOMPThreads << std::endl;
#endif
    std::cout << "Gravity: " << m_p.gravity << " m/s^2" << std::endl;
    std::cout << "Density: " << m_p.fluid.rho << " kg/m^3" << std::endl;
    std::cout << "Viscosity: " << m_p.fluid.mu << " Pa s" << std::endl;
    std::cout << "Picard relative tolerance: " << m_p.picard.relTol  << std::endl;
    std::cout << "Maximum picard iteration number: " << m_p.picard.maxIter << std::endl;
    std::cout << "End simulation time: " << m_p.time.simuTime << " s" << std::endl;
    std::cout << "Write solution every: " << m_p.time.simuDTToWrite << " s" << std::endl;
    std::cout << "Adapt time step: " << (m_p.time.adaptDT ? "yes" : "no") << std::endl;
    if(m_p.time.adaptDT)
    {
        std::cout << "Maximum time step: " << m_p.time.maxDT << " s" << std::endl;
        std::cout << "Time step reduction coefficient: " << m_p.time.coeffDTdecrease << std::endl;
        std::cout << "Time step increase coefficient: " << m_p.time.coeffDTincrease << std::endl;
    }
    else
        std::cout << "Time step: " << m_p.time.maxDT << " s" << std::endl;
}

void Solver::setInitialCondition()
{
    assert(m_mesh.getNodesNumber() != 0);

    m_mesh.setStatesNumber(3);

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        if(!m_mesh.isNodeBound(n) || m_mesh.isNodeFluidInput(n))
        {
            m_mesh.setNodeState(n, 0, m_p.initialCondition[0]);
            m_mesh.setNodeState(n, 1, m_p.initialCondition[1]);
            m_mesh.setNodeState(n, 2, m_p.initialCondition[2]);
        }
        else
        {
            m_mesh.setNodeState(n, 0, 0);
            m_mesh.setNodeState(n, 1, 0);
            m_mesh.setNodeState(n, 2, 0);
        }
    }
}

void Solver::setNodesStatesfromQ(const Eigen::VectorXd& q, unsigned short beginState, unsigned short endState)
{
    assert (q.rows() == (endState - beginState + 1)*m_mesh.getNodesNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        for (unsigned short s = beginState ; s <= endState ; ++s)
            m_mesh.setNodeState(n, s, q((s - beginState)*m_mesh.getNodesNumber() + n));
    }
}

void Solver::solveProblem()
{
    std::cout   << "================================================================"
                << std::endl
                << "                     EXECUTING THE SOLVER                       "
                << std::endl
                << "================================================================"
                << std::endl;

    displaySolverParams();

    std::cout << "----------------------------------------------------------------" << std::endl;

    setInitialCondition();

    writeData();

    m_qprev = getQFromNodesStates(0, 2);

    while(m_p.time.currentTime < m_p.time.simuTime)
    {
        if(m_verboseOutput)
        {
            std::cout << "----------------------------------------------------------------" << std::endl;
            std::cout << "Solving time step: " << m_p.time.currentTime + m_p.time.currentDT
                      << "/" << m_p.time.simuTime << " s, dt = " << m_p.time.currentDT << std::endl;
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
        else
        {
            std::cout << std::fixed << std::setprecision(3);
            std::cout << "\r" << "Solving time step: " << m_p.time.currentTime + m_p.time.currentDT
                      << "/" << m_p.time.simuTime << " s, dt = ";
            std::cout << std::scientific;
            std::cout << m_p.time.currentDT << " s" << std::flush;
        }

        if(solveCurrentTimeStep())
        {
            if(m_p.time.currentTime >= m_p.time.nextWriteTrigger)
            {
                writeData();
                m_p.time.nextWriteTrigger += m_p.time.simuDTToWrite;
            }

            if(m_p.time.adaptDT)
            {
                if(m_p.picard.currentNumIter < m_p.picard.maxIter)
                {
                    m_p.time.currentDT = std::min(m_p.time.maxDT,
                                                     m_p.time.currentDT*m_p.time.coeffDTincrease);
                }
            }

            //Remeshing step
            m_mesh.remesh();

            //We have to compute qPrev here due to new nodes !
            m_qprev = getQFromNodesStates(0, 2);
        }
        else
        {
            if(m_p.time.adaptDT)
            {
                if(m_verboseOutput)
                {
                    std::cout << "* The algorithm did not converge - Reducing the time step"
                              << std::endl;
                }

                m_p.time.currentDT /= m_p.time.coeffDTdecrease;
            }
            else
            {
                throw std::runtime_error("The solver does not seem to converge");
            }
        }

        m_p.time.currentStep++;
    }

    std::cout << std::endl;
}

bool Solver::solveCurrentTimeStep()
{
    m_mesh.saveNodesList();

    m_p.picard.currentNumIter = 0;
    Eigen::VectorXd qIter(3*m_mesh.getNodesNumber()); qIter.setZero();
    Eigen::VectorXd qIterPrev(3*m_mesh.getNodesNumber()); qIterPrev.setZero();
    double res = std::numeric_limits<double>::max();

    while(res > m_p.picard.relTol)
    {
        if(m_verboseOutput)
        {
            std::cout << "Picard algorithm (mesh position) - iteration ("
                      << m_p.picard.currentNumIter << ")" << std::endl;
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

        setNodesStatesfromQ(qIter, 0, 2);
        Eigen::VectorXd deltaPos = qIter*m_p.time.currentDT;
        m_mesh.updateNodesPositionFromSave(std::vector<double> (deltaPos.data(), deltaPos.data() + 2*m_mesh.getNodesNumber()));
        m_p.picard.currentNumIter++;

        if(m_p.picard.currentNumIter == 1)
            res = std::numeric_limits<double>::max();
        else
        {
            double num{0};
            double den{0};

            for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
            {
                if(!m_mesh.isNodeFree(n))
                {
                    num += (qIter(n) - qIterPrev(n))*(qIter(n) - qIterPrev(n));
                    den += qIterPrev(n)*qIterPrev(n);

                    num += (qIter(n + m_mesh.getNodesNumber()) - qIterPrev(n + m_mesh.getNodesNumber()))*(qIter(n + m_mesh.getNodesNumber()) - qIterPrev(n + m_mesh.getNodesNumber()));
                    den += qIterPrev(n + m_mesh.getNodesNumber())*qIterPrev(n + m_mesh.getNodesNumber());
                }
            }

            res = std::sqrt(num/den);
        }

        if(m_verboseOutput)
        {
            std::cout << " * Relative 2-norm of q: " << res << " vs "
                      << m_p.picard.relTol << std::endl;

            Eigen::VectorXd mom = m_M*(qIter.head(2*m_mesh.getNodesNumber()) - m_qprev.head(2*m_mesh.getNodesNumber()))
                                + m_K*qIter.head(2*m_mesh.getNodesNumber())
                                - Eigen::MatrixXd(m_D).transpose()*qIter.tail(m_mesh.getNodesNumber())
                                - m_F;

            Eigen::VectorXd cont = m_D*qIter.head(2*m_mesh.getNodesNumber());

            for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
            {
                if(m_mesh.isNodeBound(n) || m_mesh.isNodeFree(n))
                {
                    mom(n) = 0;
                    mom(n + m_mesh.getNodesNumber()) = 0;
                    cont(n) = 0;
                }
            }

            std::cout << " * Error on momentum: " << mom.norm() <<std::endl;
            std::cout << " * Error on mass: " << cont.norm() <<std::endl;
        }

        if(m_p.picard.currentNumIter >= m_p.picard.maxIter || std::isnan(res))
        {
            m_mesh.restoreNodesList();
            return false;
        }
    }

    m_p.time.currentTime += m_p.time.currentDT;

    return true;
}

void Solver::buildPicardSystem()
{
    assert(m_tauPSPG.size() == m_mesh.getElementNumber());

    /*A = [(matrices.M)/p.dt + matrices.K, -transpose(matrices.D);...
         (matrices.C)/p.dt + matrices.D, matrices.L];

      b = [matrices.F + (matrices.M*qPrev(1:2*Nn))/p.dt;...
           matrices.H + (matrices.C*qPrev(1:2*Nn))/p.dt];
    */

    m_A.resize(3*m_mesh.getNodesNumber(), 3*m_mesh.getNodesNumber());
    m_A.data().squeeze();
    std::vector<Eigen::Triplet<double>> indexA;
    indexA.resize(99*m_mesh.getNodesNumber());

    m_b.resize(3*m_mesh.getNodesNumber());
    m_b.setZero();

    m_M.resize(2*m_mesh.getNodesNumber(), 2*m_mesh.getNodesNumber());
    m_M.data().squeeze();
    std::vector<Eigen::Triplet<double>> indexM;
    indexM.reserve(2*3*3*m_mesh.getElementNumber());

    std::vector<Eigen::Triplet<double>> indexK;
    std::vector<Eigen::Triplet<double>> indexD;
    if(m_verboseOutput)
    {
        m_K.resize(2*m_mesh.getNodesNumber(), 2*m_mesh.getNodesNumber());
        m_K.data().squeeze();
        indexK.reserve(6*6*m_mesh.getElementNumber());

        m_D.resize(m_mesh.getNodesNumber(), 2*m_mesh.getNodesNumber());
        m_D.data().squeeze();
        indexD.reserve(3*6*m_mesh.getElementNumber());
    }

    Eigen::SparseMatrix<double> C(m_mesh.getNodesNumber(), 2*m_mesh.getNodesNumber());
    std::vector<Eigen::Triplet<double>> indexC;
    indexC.reserve(3*6*m_mesh.getElementNumber());

    m_F.resize(2*m_mesh.getNodesNumber());
    m_F.setZero();

    Eigen::VectorXd H(m_mesh.getNodesNumber());
    H.setZero();

    Eigen::MatrixXd MPrev = m_sumNTN*0.5*(m_p.fluid.rho/m_p.time.currentDT);

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        Eigen::MatrixXd Be = getB(elm);
        Eigen::MatrixXd Bep(2,3);
        Bep << Be.topLeftCorner<1,3>(),
               Be.bottomRightCorner<1,3>();

        //Me = S Nv^T Nv dV
        Eigen::MatrixXd Me = MPrev*m_mesh.getDetJ(elm);

        //Ke = S Bv^T ddev Bv dV
        Eigen::MatrixXd Ke = 0.5*Be.transpose()*m_ddev*Be*m_mesh.getDetJ(elm); //same matrices for all gauss point ^^

        //De = S Np^T m Bv dV
        Eigen::MatrixXd De(3,6); De.setZero();

        for (unsigned short k = 0 ; k < m_N.size() ; ++k)
            De += (m_N[k].topLeftCorner<1,3>()).transpose()*m_m.transpose()*Be*GP2Dweight<double>[k];

        De *= 0.5*m_mesh.getDetJ(elm);

        //Ce = tauPSPG S Bp^T Nv dV
        Eigen::MatrixXd Ce(3,6); Ce.setZero();

        for (unsigned short k = 0 ; k < m_N.size() ; ++k)
            Ce += Bep.transpose()*m_N[k]*GP2Dweight<double>[k];

        Ce *= (0.5*m_tauPSPG[elm]/m_p.time.currentDT)*m_mesh.getDetJ(elm);

        //Le = tauPSPG S Bp^T Bp dV
        Eigen::MatrixXd Le = 0.5*(m_tauPSPG[elm]/m_p.fluid.rho)
                             *Bep.transpose()*Bep
                             *m_mesh.getDetJ(elm);

        const Eigen::Vector2d b(0, -m_p.gravity);

        //Fe = S Nv^T bodyforce dV
        Eigen::MatrixXd Fe(6,1); Fe.setZero();

        for (unsigned short k = 0 ; k < m_N.size() ; ++k)
            Fe += m_N[k].transpose()*b*GP2Dweight<double>[k];

        Fe *= 0.5*m_p.fluid.rho*m_mesh.getDetJ(elm);

        //He = tauPSPG S Bp^T bodyforce dV
        Eigen::VectorXd He = 0.5*m_tauPSPG[elm]*Bep.transpose()
                             *b*m_mesh.getDetJ(elm);

        //Big push_back ^^
        #pragma omp critical
        {
            for(unsigned short i = 0 ; i < 3 ; ++i)
            {
                for(unsigned short j = 0 ; j < 3 ; ++j)
                {
                    /********************************************************************
                                                 Build M/dt
                    ********************************************************************/
                    indexM.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j],
                                                            Me(i,j)));

                    indexM.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j] + m_mesh.getNodesNumber(),
                                                            Me(i+3,j+3)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j],
                                                            Me(i,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j] + m_mesh.getNodesNumber(),
                                                            Me(i+3,j+3)));

                    /********************************************************************
                                                    Build K
                    ********************************************************************/
                    if(m_verboseOutput)
                    {
                        indexK.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                                m_mesh.getElement(elm)[j],
                                                                Ke(i,j)));

                        indexK.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                                m_mesh.getElement(elm)[j] + m_mesh.getNodesNumber(),
                                                                Ke(i,j+3)));

                        indexK.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.getNodesNumber(),
                                                                m_mesh.getElement(elm)[j],
                                                                Ke(i+3,j)));

                        indexK.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.getNodesNumber(),
                                                                m_mesh.getElement(elm)[j] + m_mesh.getNodesNumber(),
                                                                Ke(i+3,j+3)));
                    }

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j],
                                                            Ke(i,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j] + m_mesh.getNodesNumber(),
                                                            Ke(i,j+3)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j],
                                                            Ke(i+3,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j] + m_mesh.getNodesNumber(),
                                                            Ke(i+3,j+3)));

                    /********************************************************************
                                                    Build D
                    ********************************************************************/
                    if(m_verboseOutput)
                    {
                        indexD.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                                m_mesh.getElement(elm)[j],
                                                                De(i,j)));

                        indexD.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                                m_mesh.getElement(elm)[j] + m_mesh.getNodesNumber(),
                                                                De(i,j+3)));
                    }

                    //Part D of A
                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j],
                                                            De(i,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j] + m_mesh.getNodesNumber(),
                                                            De(i,j+3)));

                    //Part -D^T of A
                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[j],
                                                            m_mesh.getElement(elm)[i] + 2*m_mesh.getNodesNumber(),
                                                            -De(i,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[j] + m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[i] + 2*m_mesh.getNodesNumber(),
                                                            -De(i,j+3)));

                    /********************************************************************
                                                Build C/dt
                    ********************************************************************/
                    indexC.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j],
                                                            Ce(i,j)));

                    indexC.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j] + m_mesh.getNodesNumber(),
                                                            Ce(i,j+3)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j],
                                                            Ce(i,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j] + m_mesh.getNodesNumber(),
                                                            Ce(i,j+3)));

                    /********************************************************************
                                                Build L
                    ********************************************************************/
                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.getNodesNumber(),
                                                            m_mesh.getElement(elm)[j] + 2*m_mesh.getNodesNumber(),
                                                            Le(i,j)));
                }

                /************************************************************************
                                                Build f
                ************************************************************************/
                m_F(m_mesh.getElement(elm)[i]) += Fe(i);

                m_F(m_mesh.getElement(elm)[i] + m_mesh.getNodesNumber()) +=
                Fe(i + 3);


                /************************************************************************
                                                Build f
                ************************************************************************/
                H(m_mesh.getElement(elm)[i]) += He(i);
            }
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

    m_b << m_F + m_M*m_qprev.head(2*m_mesh.getNodesNumber()), H + C*m_qprev.head(2*m_mesh.getNodesNumber());
}

void Solver::computeTauPSPG()
{
    m_tauPSPG.resize(m_mesh.getElementNumber());

    double U = 0;
    std::size_t trueNnodes = 0;
    #pragma omp parallel for default(shared) reduction(+: U, trueNnodes)
    for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        if(!m_mesh.isNodeFree(n))
        {
            U += std::sqrt(m_mesh.getNodeState(n, 0)*m_mesh.getNodeState(n, 0) +
                           m_mesh.getNodeState(n, 1)*m_mesh.getNodeState(n, 1));

            trueNnodes++;
        }
    }
    U /= static_cast<double>(trueNnodes);

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        const double h = std::sqrt(2*m_mesh.getDetJ(elm)/M_PI);

        m_tauPSPG[elm] = 1/std::sqrt((2/m_p.time.currentDT)*(2/m_p.time.currentDT)
                                 + (2*U/h)*(2*U/h)
                                 + 9*(4*m_p.fluid.mu/(h*h*m_p.fluid.rho))*(4*m_p.fluid.mu/(h*h*m_p.fluid.rho)));
    }
}

void Solver::writeData() const
{
    gmsh::model::add("theModel");
    gmsh::model::setCurrent("theModel");
    if(m_p.writeAs == "Nodes" || m_p.writeAs == "NodesElements")
        gmsh::model::addDiscreteEntity(0, 1);
    if(m_p.writeAs == "Elements" || m_p.writeAs == "NodesElements")
        gmsh::model::addDiscreteEntity(2, 2);

    if(m_p.whatToWrite[0])
        gmsh::view::add("u", 1);

    if(m_p.whatToWrite[1])
        gmsh::view::add("v", 2);

    if(m_p.whatToWrite[2])
        gmsh::view::add("p", 3);

    if(m_p.whatToWrite[3])
        gmsh::view::add("ke", 4);

    if(m_p.whatToWrite[4])
        gmsh::view::add("velocity", 5);

    std::vector<std::size_t> nodesTags(m_mesh.getNodesNumber());
    std::vector<double> nodesCoord(3*m_mesh.getNodesNumber());
    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        nodesTags[n] = n + 1;
        nodesCoord[3*n] = m_mesh.getNodePosition(n, 0);
        nodesCoord[3*n + 1] = m_mesh.getNodePosition(n, 1);
        nodesCoord[3*n + 2] = 0;
    }

    gmsh::model::mesh::addNodes(0, 1, nodesTags, nodesCoord);

    if(m_p.writeAs == "Elements" || m_p.writeAs == "NodesElements")
    {
        std::vector<std::size_t> elementTags(m_mesh.getElementNumber());
        std::vector<std::size_t> nodesTagsPerElement(3*m_mesh.getElementNumber());
        #pragma omp parallel for default(shared)
        for(std::size_t i = 0 ; i < m_mesh.getElementNumber() ; ++i)
        {
            elementTags[i] = i + 1;
            nodesTagsPerElement[3*i] = m_mesh.getElement(i)[0] + 1;
            nodesTagsPerElement[3*i + 1] = m_mesh.getElement(i)[1] + 1;
            nodesTagsPerElement[3*i + 2] = m_mesh.getElement(i)[2] + 1;
        }

        gmsh::model::mesh::addElementsByType(2, 2, elementTags, nodesTagsPerElement);
    }

    if(m_p.writeAs == "Nodes" || m_p.writeAs == "NodesElements")
        gmsh::model::mesh::addElementsByType(1, 15, nodesTags, nodesTags);

    std::vector<std::vector<double>> dataU;
    std::vector<std::vector<double>> dataV;
    std::vector<std::vector<double>> dataP;
    std::vector<std::vector<double>> dataKe;
    std::vector<std::vector<double>> dataVelocity;

    if(m_p.whatToWrite[0])
        dataU.resize(m_mesh.getNodesNumber());
    if(m_p.whatToWrite[1])
        dataV.resize(m_mesh.getNodesNumber());
    if(m_p.whatToWrite[2])
        dataP.resize(m_mesh.getNodesNumber());
    if(m_p.whatToWrite[3])
        dataKe.resize(m_mesh.getNodesNumber());
    if(m_p.whatToWrite[4])
        dataVelocity.resize(m_mesh.getNodesNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        if(m_p.whatToWrite[0])
        {
            const std::vector<double> u{m_mesh.getNodeState(n, 0)};
            dataU[n] = u;
        }
        if(m_p.whatToWrite[1])
        {
            const std::vector<double> v{m_mesh.getNodeState(n, 1)};
            dataV[n] = v;
        }
        if(m_p.whatToWrite[2])
        {
            const std::vector<double> p{m_mesh.getNodeState(n, 2)};
            dataP[n] = p;
        }
        if(m_p.whatToWrite[3])
        {
            const std::vector<double> ke{0.5*(m_mesh.getNodeState(n, 0)*m_mesh.getNodeState(n, 0)
                                            + m_mesh.getNodeState(n, 1)*m_mesh.getNodeState(n, 1))};
            dataKe[n] = ke;
        }
        if(m_p.whatToWrite[4])
        {
            const std::vector<double> velocity{m_mesh.getNodeState(n, 0), m_mesh.getNodeState(n, 1), 0};;
            dataVelocity[n] = velocity;
        }
    }

    if(m_p.whatToWrite[0])
        gmsh::view::addModelData(1, m_p.time.currentStep, "theModel", "NodeData", nodesTags, dataU, m_p.time.currentTime, 1);
    if(m_p.whatToWrite[1])
        gmsh::view::addModelData(2, m_p.time.currentStep, "theModel", "NodeData", nodesTags, dataV, m_p.time.currentTime, 1);
    if(m_p.whatToWrite[2])
        gmsh::view::addModelData(3, m_p.time.currentStep, "theModel", "NodeData", nodesTags, dataP, m_p.time.currentTime, 1);
    if(m_p.whatToWrite[3])
        gmsh::view::addModelData(4, m_p.time.currentStep, "theModel", "NodeData", nodesTags, dataKe, m_p.time.currentTime, 1);
    if(m_p.whatToWrite[4])
        gmsh::view::addModelData(5, m_p.time.currentStep, "theModel", "NodeData", nodesTags, dataVelocity, m_p.time.currentTime, 3);

    const std::string baseName = m_p.resultsName.substr(0, m_p.resultsName.find(".msh"));

    for(unsigned short i = 0; i < m_p.whatToWrite.size(); ++i)
    {
        if(m_p.whatToWrite[i] == true)
            gmsh::view::write(i + 1, baseName + "_" + std::to_string(m_p.time.currentTime) + ".msh" , true);
    }

    gmsh::model::remove();
}
