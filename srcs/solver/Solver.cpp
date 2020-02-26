#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include <gmsh.h>
#include <nlohmann/json.hpp>

#include "Solver.hpp"
#include "../quadrature/gausslegendre.hpp"

Solver::Solver(const Params& params, Mesh& mesh, std::string resultsName) :
m_mesh(mesh)
{
    nlohmann::json j = params.getJSON();

    m_verboseOutput = j["verboseOutput"].get<bool>();

    m_p.hchar = j["RemeshingParams"]["hchar"].get<double>();

    m_p.gravity = j["SolverIncompressibleParams"]["gravity"].get<double>();

    m_p.fluid.rho = j["SolverIncompressibleParams"]["FluidParams"]["rho"].get<double>();
    m_p.fluid.mu = j["SolverIncompressibleParams"]["FluidParams"]["mu"].get<double>();

    m_p.picard.relTol = j["SolverIncompressibleParams"]["PicardParams"]["relTol"].get<double>();
    m_p.picard.maxIter = j["SolverIncompressibleParams"]["PicardParams"]["maxIter"].get<unsigned int>();

    m_p.time.adaptDT = j["SolverIncompressibleParams"]["TimeParams"]["adaptDT"].get<bool>();
    m_p.time.coeffDTincrease = j["SolverIncompressibleParams"]["TimeParams"]["coeffDTincrease"].get<double>();
    m_p.time.coeffDTdecrease = j["SolverIncompressibleParams"]["TimeParams"]["coeffDTdecrease"].get<double>();
    m_p.time.maxDT = j["SolverIncompressibleParams"]["TimeParams"]["maxDT"].get<double>();
    m_p.time.simuTime = j["SolverIncompressibleParams"]["TimeParams"]["simuTime"].get<double>();
    m_p.time.simuDTToWrite = j["SolverIncompressibleParams"]["TimeParams"]["simuDTToWrite"].get<double>();
    m_p.time.currentDT = m_p.time.maxDT;
    m_p.time.nextWriteTrigger = m_p.time.simuDTToWrite;

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

    gmsh::initialize();

#ifndef NDEBUG
    gmsh::option::setNumber("General.Terminal", 1);
#else
    gmsh::option::setNumber("General.Terminal", 0);
#endif // DEBUG

    gmsh::option::setNumber("General.NumThreads", params.getNumOMPThreads());

    m_p.resultsName = resultsName;
}

Solver::~Solver()
{
    gmsh::finalize();
}

void Solver::solveProblem()
{
    std::cout   << "================================================================"
                << std::endl
                << "                     EXECUTING THE SOLVER                       "
                << std::endl
                << "================================================================"
                << std::endl;

    std::cout << "Eigen sparse solver: SparseLU" << std::endl;
    std::cout << "Gravity: " << m_p.gravity << " m/s^2" << std::endl;
    std::cout << "Density: " << m_p.fluid.rho << " kg/m^3" << std::endl;
    std::cout << "Viscosity: " << m_p.fluid.mu << " Pa s" << std::endl;
    std::cout << "Picard relative tolerance: " << m_p.picard.relTol  << std::endl;
    std::cout << "Maxixum picard iteration number: " << m_p.picard.maxIter << std::endl;
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

    std::cout << "----------------------------------------------------------------" << std::endl;

    if(!m_verboseOutput)
        std::cout << std::fixed << std::setprecision(3);

    double time = 0;

    unsigned int step = 0;

    m_mesh.remesh();

    m_qprev.resize(3*m_mesh.nodesList.size());
    for(std::size_t n = 0 ; n < m_mesh.nodesList.size() ; ++n)
    {
        m_qprev(n) = m_mesh.nodesList[n].states[0];
        m_qprev(n + m_mesh.nodesList.size()) = m_mesh.nodesList[n].states[1];
        m_qprev(n + 2*m_mesh.nodesList.size()) = m_mesh.nodesList[n].states[2];
    }

    _write(time, step);

    while(time < m_p.time.simuTime)
    {
        if(m_verboseOutput)
        {
            std::cout << "----------------------------------------------------------------" << std::endl;
            std::cout << "Solving time step: " << time + m_p.time.currentDT
                      << "/" << m_p.time.simuTime << " s" << std::endl;
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
        else
        {
            std::cout << "\r" << "Solving time step: " << time + m_p.time.currentDT
                      << "/" << m_p.time.simuTime << " s" << std::flush;
        }

        if(_solveSystem())
        {
            time = time + m_p.time.currentDT;

            #pragma omp parallel for default(none) shared(m_mesh, m_q, m_p)
            for(std::size_t n = 0 ; n < m_mesh.nodesList.size() ; ++n)
            {
                m_mesh.nodesList[n].states[0] = m_q(n);
                m_mesh.nodesList[n].states[1] = m_q(n + m_mesh.nodesList.size());
                m_mesh.nodesList[n].states[2] = m_q(n + 2*m_mesh.nodesList.size());

                m_mesh.nodesList[n].position[0] += m_mesh.nodesList[n].states[0]*m_p.time.currentDT;
                m_mesh.nodesList[n].position[1] += m_mesh.nodesList[n].states[1]*m_p.time.currentDT;
            }

            if(time >= m_p.time.nextWriteTrigger)
            {
                _write(time, step);
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

            if(time + m_p.time.currentDT > m_p.time.simuTime)
                        m_p.time.currentDT = m_p.time.simuTime - time;

            if(m_mesh.removeNodes())
                m_mesh.remesh();

            m_mesh.addNodes();
            m_mesh.remesh();

            m_qprev.resize(3*m_mesh.nodesList.size());

            #pragma omp parallel for default(none) shared(m_mesh, m_qprev)
            for(std::size_t n = 0 ; n < m_mesh.nodesList.size() ; ++n)
            {
                m_qprev(n) = m_mesh.nodesList[n].states[0];
                m_qprev(n + m_mesh.nodesList.size()) = m_mesh.nodesList[n].states[1];
                m_qprev(n + 2*m_mesh.nodesList.size()) = m_mesh.nodesList[n].states[2];
            }
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

        step++;
    }

    _write(time, step);

    std::cout << std::endl;

    //_writeFinalize();
}

bool Solver::_solveSystem()
{
    std::vector<Node> nodesListSave {m_mesh.nodesList};

    std::vector<bool> indices = _getIndices();

    m_p.picard.currentNumIter = 0;
    Eigen::VectorXd qIter(3*m_mesh.nodesList.size()); qIter.setZero();
    Eigen::VectorXd qIterPrev(3*m_mesh.nodesList.size()); qIterPrev.setZero();
    double res = std::numeric_limits<double>::max();

    while(res > m_p.picard.relTol)
    {
        if(m_verboseOutput)
        {
            std::cout << "Picard algorithm (mesh position) - iteration ("
                      << m_p.picard.currentNumIter << ")" << std::endl;
        }

        qIterPrev = qIter;

        m_mesh.computeDetJ();
        m_mesh.computeInvJ();
        _computeTauPSPG();
        _buildPicardSystem();

        m_A.prune([indices](int i, int j, float)
                  {
                      return !indices[i];
                  });

        for(std::size_t i = 0 ; i < indices.size() ; ++i)
        {
            if(indices[i])
            {
                m_b(i) = m_qprev(i);
                m_A.coeffRef(i, i) = 1;
            }
        }

        m_solverLU.compute(m_A);
        if(m_solverLU.info() == Eigen::Success)
            qIter = m_solverLU.solve(m_b);
        else
        {
            m_mesh.nodesList = nodesListSave;
            return false;
        }

        #pragma omp parallel for default(none) shared(m_mesh, qIter, m_p, nodesListSave)
        for(std::size_t n = 0 ; n < m_mesh.nodesList.size() ; ++n)
        {
            m_mesh.nodesList[n].states[0] = qIter[n];
            m_mesh.nodesList[n].states[1] = qIter[n + m_mesh.nodesList.size()];
            m_mesh.nodesList[n].states[2] = qIter[n + 2*m_mesh.nodesList.size()];
            m_mesh.nodesList[n].position[0] = nodesListSave[n].position[0]
                                            + m_p.time.currentDT*m_mesh.nodesList[n].states[0];
            m_mesh.nodesList[n].position[1] = nodesListSave[n].position[1]
                                            + m_p.time.currentDT*m_mesh.nodesList[n].states[1];
        }

        m_p.picard.currentNumIter++;

        if(m_p.picard.currentNumIter == 1)
            res = std::numeric_limits<double>::max();
        else
        {
            double num{0};
            double den{0};

            for(std::size_t i = 0 ; i < 2*m_mesh.nodesList.size() ; ++i)
            {
                if(!indices[i])
                {
                    num += (qIter(i) - qIterPrev(i))*(qIter(i) - qIterPrev(i));
                    den += qIterPrev(i)*qIterPrev(i);
                }
            }

            res = std::sqrt(num/den);
        }

        if(m_verboseOutput)
        {
            std::cout << " * Relative 2-norm of q: " << res << " vs "
                      << m_p.picard.relTol << std::endl;

            Eigen::VectorXd mom = m_M*(qIter.head(2*m_mesh.nodesList.size()) - m_qprev.head(2*m_mesh.nodesList.size()))
                                + m_K*qIter.head(2*m_mesh.nodesList.size())
                                - Eigen::MatrixXd(m_D).transpose()*qIter.tail(m_mesh.nodesList.size())
                                - m_F;

            Eigen::VectorXd cont = m_D*qIter.head(2*m_mesh.nodesList.size());

            for(std::size_t i = 0 ; i < 2*m_mesh.nodesList.size() ; ++i)
            {
                if(indices[i])
                {
                    mom(i) = 0;
                    if(i < m_mesh.nodesList.size())
                        cont(i) = 0;
                }
            }

            std::cout << " * Error on momentum: " << mom.norm() <<std::endl;
            std::cout << " * Error on mass: " << cont.norm() <<std::endl;
        }

        if(m_p.picard.currentNumIter >= m_p.picard.maxIter || std::isnan(res))
        {
            m_mesh.nodesList = nodesListSave;
            return false;
        }
    }

    m_q = qIter;
    m_mesh.nodesList = nodesListSave;

    #pragma omp parallel for default(none) shared(m_mesh, m_q)
    for(std::size_t n = 0 ; n < m_mesh.nodesList.size() ; ++n)
    {
        if(m_mesh.nodesList[n].isFree && !m_mesh.nodesList[n].isBound)
        {
            m_q(m_mesh.nodesList.size() + n) = m_q (m_mesh.nodesList.size() + n) -
                                                    m_p.time.currentDT*m_p.gravity*m_p.fluid.rho*
                                                    m_p.hchar*m_p.hchar*0.5;
        }
    }

    return true;
}

void Solver::_buildPicardSystem()
{
    assert(!m_tauPSPG.empty());

    /*A = [(matrices.M)/p.dt + matrices.K, -transpose(matrices.D);...
         (matrices.C)/p.dt + matrices.D, matrices.L];

      b = [matrices.F + (matrices.M*qPrev(1:2*Nn))/p.dt;...
           matrices.H + (matrices.C*qPrev(1:2*Nn))/p.dt];
    */

    std::vector<Eigen::MatrixXd> N = m_mesh.getN();
    std::vector<std::vector<Eigen::MatrixXd>> B(m_mesh.getElementNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        B[elm] = m_mesh.getB(elm);
    }

    assert(N.size() == 3);

    const Eigen::MatrixXd BrightP = (Eigen::MatrixXd(6,3) << 1, 0, 0,
                                                             0, 1, 0,
                                                             0, 0, 1,
                                                             1, 0, 0,
                                                             0, 1, 0,
                                                             0, 0, 1).finished();

    const Eigen::MatrixXd BleftP = (Eigen::MatrixXd(2,3) << 1, 0, 0,
                                                            0, 1, 0).finished();

    std::vector<Eigen::Triplet<double>> index;

    m_A.resize(3*m_mesh.nodesList.size(), 3*m_mesh.nodesList.size());
    m_A.data().squeeze();
    std::vector<Eigen::Triplet<double>> indexA;
    indexA.resize(99*m_mesh.nodesList.size());

    m_b.resize(3*m_mesh.nodesList.size());
    m_b.setZero();

    m_M.resize(2*m_mesh.nodesList.size(), 2*m_mesh.nodesList.size());
    m_M.data().squeeze();
    std::vector<Eigen::Triplet<double>> indexM;
    indexM.reserve(2*3*3*m_mesh.getElementNumber());

    std::vector<Eigen::Triplet<double>> indexK;
    std::vector<Eigen::Triplet<double>> indexD;
    if(m_verboseOutput)
    {
        m_K.resize(2*m_mesh.nodesList.size(), 2*m_mesh.nodesList.size());
        m_K.data().squeeze();
        indexK.reserve(6*6*m_mesh.getElementNumber());

        m_D.resize(m_mesh.nodesList.size(), 2*m_mesh.nodesList.size());
        m_D.data().squeeze();
        indexD.reserve(3*6*m_mesh.getElementNumber());
    }

    m_C.resize(m_mesh.nodesList.size(), 2*m_mesh.nodesList.size());
    m_C.data().squeeze();
    std::vector<Eigen::Triplet<double>> indexC;
    indexC.reserve(3*6*m_mesh.getElementNumber());

    m_F.resize(2*m_mesh.nodesList.size());
    m_F.setZero();

    m_H.resize(m_mesh.nodesList.size());
    m_H.setZero();

    Eigen::MatrixXd MPrev(6, 6); MPrev.setZero();
    for (unsigned short k = 0 ; k < N.size() ; ++k)
        MPrev += N[k].transpose()*N[k]*GP2Dweight<double>[k];

    MPrev *= 0.5*(m_p.fluid.rho/m_p.time.currentDT);

    const Eigen::DiagonalMatrix<double, 3> ddev(2, 2, 1);

    const Eigen::Vector3d m(1, 1, 0);

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        Eigen::MatrixXd Me = MPrev*m_mesh.getDetJ(elm);

        Eigen::MatrixXd Ke = 0.5*m_p.fluid.mu*B[elm][0].transpose()*ddev*B[elm][0]*m_mesh.getDetJ(elm); //same matrices for all gauss point ^^

        Eigen::MatrixXd De(3,6); De.setZero();

        for (unsigned short k = 0 ; k < N.size() ; ++k)
            De += (N[k].topLeftCorner<1,3>()).transpose()*m.transpose()*B[elm][k]*GP2Dweight<double>[k];

        De *= 0.5*m_mesh.getDetJ(elm);

        Eigen::MatrixXd Ce(3,6); Ce.setZero();

        for (unsigned short k = 0 ; k < N.size() ; ++k)
            Ce +=(BleftP*B[elm][k]*BrightP).transpose()*N[k]*GP2Dweight<double>[k];

        Ce *= (0.5*m_tauPSPG[elm]/m_p.time.currentDT)*m_mesh.getDetJ(elm);

        Eigen::MatrixXd Le = 0.5*(m_tauPSPG[elm]/m_p.fluid.rho)
                             *(BleftP*B[elm][0]*BrightP).transpose()*(BleftP*B[elm][0]*BrightP)
                             *m_mesh.getDetJ(elm);

        const Eigen::Vector2d b(0, -m_p.gravity);

        Eigen::MatrixXd Fe(6,1); Fe.setZero();

        for (unsigned short k = 0 ; k < N.size() ; ++k)
            Fe += N[k].transpose()*b*GP2Dweight<double>[k];

        Fe *= 0.5*m_p.fluid.rho*m_mesh.getDetJ(elm);

        Eigen::VectorXd He = 0.5*m_tauPSPG[elm]*(BleftP*B[elm][0]*BrightP).transpose()
                             *b*m_mesh.getDetJ(elm);

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

                    indexM.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                            m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                            Me(i+3,j+3)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j],
                                                            Me(i,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                            m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
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
                                                                m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                                Ke(i,j+3)));

                        indexK.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                                m_mesh.getElement(elm)[j],
                                                                Ke(i+3,j)));

                        indexK.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                                m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                                Ke(i+3,j+3)));
                    }

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j],
                                                            Ke(i,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                            Ke(i,j+3)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                            m_mesh.getElement(elm)[j],
                                                            Ke(i+3,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                            m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
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
                                                                m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                                De(i,j+3)));
                    }

                    //Part D of A
                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                            m_mesh.getElement(elm)[j],
                                                            De(i,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                            m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                            De(i,j+3)));

                    //Part -D^T of A
                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[j],
                                                            m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                            -De(i,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                            m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                            -De(i,j+3)));

                    /********************************************************************
                                                Build C/dt
                    ********************************************************************/
                    indexC.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j],
                                                            Ce(i,j)));

                    indexC.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                            m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                            Ce(i,j+3)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                            m_mesh.getElement(elm)[j],
                                                            Ce(i,j)));

                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                            m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                            Ce(i,j+3)));

                    /********************************************************************
                                                Build L
                    ********************************************************************/
                    indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                            m_mesh.getElement(elm)[j] + 2*m_mesh.nodesList.size(),
                                                            Le(i,j)));
                }

                /************************************************************************
                                                Build f
                ************************************************************************/
                m_F(m_mesh.getElement(elm)[i]) += Fe(i);

                m_F(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size()) +=
                Fe(i + 3);


                /************************************************************************
                                                Build f
                ************************************************************************/
                m_H(m_mesh.getElement(elm)[i]) += He(i);

            }
        }
    }

    m_M.setFromTriplets(indexM.begin(), indexM.end());

    if(m_verboseOutput)
    {
        m_K.setFromTriplets(indexK.begin(), indexK.end());
        m_D.setFromTriplets(indexD.begin(), indexD.end());
    }

    m_C.setFromTriplets(indexC.begin(), indexC.end());

    /********************************************************************************
                                        Compute A and b
    ********************************************************************************/
    m_A.setFromTriplets(indexA.begin(), indexA.end());

    m_b << m_F + m_M*m_qprev.head(2*m_mesh.nodesList.size()), m_H + m_C*m_qprev.head(2*m_mesh.nodesList.size());
}

void Solver::_computeTauPSPG()
{
    m_tauPSPG.resize(m_mesh.getElementNumber());

    double U = 0;
    std::size_t trueNnodes = 0;
    #pragma omp parallel for default(shared) reduction(+: U, trueNnodes)
    for(std::size_t n = 0 ; n < m_mesh.nodesList.size() ; ++n)
    {
        if(m_mesh.nodesList[n].isBound == false ||
           (m_mesh.nodesList[n].isBound == true && m_mesh.nodesList[n].isFree == false))
        {
            U += std::sqrt(m_mesh.nodesList[n].states[0]*
               m_mesh.nodesList[n].states[0] +
               m_mesh.nodesList[n].states[1]*
               m_mesh.nodesList[n].states[1]);

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

std::vector<bool> Solver::_getIndices() const
{
    std::vector<bool> indices(3*m_mesh.nodesList.size(), false);

    #pragma omp parallel for default(shared)
    for (std::size_t i = 0 ; i < m_mesh.nodesList.size() ; ++i)
    {
        if(m_mesh.nodesList[i].isBound || m_mesh.nodesList[i].isFree)
        {
            indices[i] = true;
            indices[i + m_mesh.nodesList.size()] = true;
        }

        if(m_mesh.nodesList[i].isFree)
            indices[i + 2*m_mesh.nodesList.size()] = true;
    }

    return indices;
}

void Solver::_write(double time, unsigned int step)
{
    gmsh::model::add("theModel");
    gmsh::model::setCurrent("theModel");
    gmsh::model::addDiscreteEntity(m_mesh.getMeshDim(), 1);

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

    std::vector<std::size_t> nodesTags(m_mesh.nodesList.size());
    std::vector<double> nodesCoord(3*m_mesh.nodesList.size());
    #pragma omp parallel for default(shared)
    for(std::size_t i = 0 ; i < m_mesh.nodesList.size() ; ++i)
    {
        nodesTags[i] = i + 1;
        nodesCoord[3*i] = m_mesh.nodesList[i].position[0];
        nodesCoord[3*i + 1] = m_mesh.nodesList[i].position[1];
        nodesCoord[3*i + 2] = 0;
    }

    gmsh::model::mesh::addNodes(m_mesh.getMeshDim(), 1, nodesTags, nodesCoord);

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

    gmsh::model::mesh::addElementsByType(1, 2, elementTags, nodesTagsPerElement);

    std::vector<std::vector<double>> dataU;
    std::vector<std::vector<double>> dataV;
    std::vector<std::vector<double>> dataP;
    std::vector<std::vector<double>> dataKe;
    std::vector<std::vector<double>> dataVelocity;

    if(m_p.whatToWrite[0])
        dataU.resize(m_mesh.nodesList.size());
    if(m_p.whatToWrite[1])
        dataV.resize(m_mesh.nodesList.size());
    if(m_p.whatToWrite[2])
        dataP.resize(m_mesh.nodesList.size());
    if(m_p.whatToWrite[3])
        dataKe.resize(m_mesh.nodesList.size());
    if(m_p.whatToWrite[4])
        dataVelocity.resize(m_mesh.nodesList.size());

    #pragma omp parallel for default(shared)
    for(std::size_t i = 0 ; i < m_mesh.nodesList.size() ; ++i)
    {
        if(m_p.whatToWrite[0])
        {
            const std::vector<double> u{m_mesh.nodesList[i].states[0]};
            dataU[i] = u;
        }
        if(m_p.whatToWrite[1])
        {
            const std::vector<double> v{m_mesh.nodesList[i].states[1]};
            dataV[i] = v;
        }
        if(m_p.whatToWrite[2])
        {
            const std::vector<double> p{m_mesh.nodesList[i].states[2]};
            dataP[i] = p;
        }
        if(m_p.whatToWrite[3])
        {
            const std::vector<double> ke{0.5*(m_mesh.nodesList[i].states[0]*m_mesh.nodesList[i].states[0]
                                        +m_mesh.nodesList[i].states[1]*m_mesh.nodesList[i].states[1])};
            dataKe[i] = ke;
        }
        if(m_p.whatToWrite[4])
        {
            const std::vector<double> velocity{m_mesh.nodesList[i].states[0], m_mesh.nodesList[i].states[1], 0};;
            dataVelocity[i] = velocity;
        }
    }

    if(m_p.whatToWrite[0])
        gmsh::view::addModelData(1, step, "theModel", "NodeData", nodesTags, dataU, time, 1);
    if(m_p.whatToWrite[1])
        gmsh::view::addModelData(2, step, "theModel", "NodeData", nodesTags, dataV, time, 1);
    if(m_p.whatToWrite[2])
        gmsh::view::addModelData(3, step, "theModel", "NodeData", nodesTags, dataP, time, 1);
    if(m_p.whatToWrite[3])
        gmsh::view::addModelData(4, step, "theModel", "NodeData", nodesTags, dataKe, time, 1);
    if(m_p.whatToWrite[4])
        gmsh::view::addModelData(5, step, "theModel", "NodeData", nodesTags, dataVelocity, time, 3);

    const std::string baseName = m_p.resultsName.substr(0, m_p.resultsName.find(".msh"));

    for(unsigned short i = 0; i < m_p.whatToWrite.size(); ++i)
    {
        if(m_p.whatToWrite[i] == true)
            gmsh::view::write(i + 1, baseName + "_" + std::to_string(time) + ".msh" , true);
    }

    gmsh::model::remove();
}

void Solver::_writeFinalize()
{
}
