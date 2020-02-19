#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>

#include "Solver.hpp"
#include "../quadrature/gausslegendre.hpp"

Solver::Solver(const Params& params, Mesh& mesh, ProblemType problemType) :
m_params(params), m_mesh(mesh), m_currentDT(params.maxTimeStep)
{
    if(problemType == INCOMPRESSIBLE_PSPG)
    {
        m_problemType = problemType;
    }
    else
        throw std::runtime_error("Unsupported problem type.");
}

Solver::~Solver()
{

}

void Solver::solveProblem()
{
    double time = 0;

    m_mesh.remesh();

    m_qprev.resize(3*m_mesh.nodesList.size());
    for(std::size_t n = 0 ; n < m_mesh.nodesList.size() ; ++n)
    {
        m_qprev(n) = m_mesh.nodesList[n].states[0];
        m_qprev(n + m_mesh.nodesList.size()) = m_mesh.nodesList[n].states[1];
        m_qprev(n + 2*m_mesh.nodesList.size()) = m_mesh.nodesList[n].states[2];
    }

    while(time < m_params.simuTime)
    {
        std::cout << "================================================" << std::endl;
        std::cout << "Solving time step: " << time << "/" << m_params.simuTime << std::endl;
        std::cout << "================================================" << std::endl;

        if(_solveSystem())
        {
            time = time + m_currentDT;
            for(std::size_t n = 0 ; n < m_mesh.nodesList.size() ; ++n)
            {
                m_mesh.nodesList[n].states[0] = m_q(n);
                m_mesh.nodesList[n].states[1] = m_q(n + m_mesh.nodesList.size());
                m_mesh.nodesList[n].states[2] = m_q(n + 2*m_mesh.nodesList.size());

                m_mesh.nodesList[n].position[0] += m_mesh.nodesList[n].states[0]*m_currentDT;
                m_mesh.nodesList[n].position[1] += m_mesh.nodesList[n].states[1]*m_currentDT;
            }

            if(m_params.adaptDT)
            {
                if(m_numIterSolverMsh < m_params.sideParams[1])
                {
                    m_currentDT = std::min(m_params.maxTimeStep,
                                                     m_currentDT*m_params.sideParams[2]);
                }
            }

            if(m_mesh.removeNodes())
                m_mesh.remesh();

            m_mesh.addNodes();
            m_mesh.remesh();

            m_qprev.resize(3*m_mesh.nodesList.size());
            for(std::size_t n = 0 ; n < m_mesh.nodesList.size() ; ++n)
            {
                m_qprev(n) = m_mesh.nodesList[n].states[0];
                m_qprev(n + m_mesh.nodesList.size()) = m_mesh.nodesList[n].states[1];
                m_qprev(n + 2*m_mesh.nodesList.size()) = m_mesh.nodesList[n].states[2];
            }
        }
        else
        {
            if(m_params.adaptDT)
            {
                std::cout << "The algorithm did not converge - Reducing the time step"
                          << std::endl;

                m_currentDT = m_currentDT/m_params.sideParams[3];
            }
            else
            {
                throw std::runtime_error("The solver does not seem to converge");
            }
        }
    }
}

bool Solver::_solveSystem()
{
    std::vector<Node> nodesListSave {m_mesh.nodesList};

    std::vector<bool> indices = _getIndices();

    m_numIterSolverMsh = 0;
    Eigen::VectorXd qIter(3*m_mesh.nodesList.size()); qIter.setZero();
    Eigen::VectorXd qIterPrev(3*m_mesh.nodesList.size()); qIterPrev.setZero();
    double res = std::numeric_limits<double>::max();

    while(res > m_params.sideParams[0])
    {
        std::cout << "Picard algorithm (mesh position) - iteration ("
                  << m_numIterSolverMsh << ")" << std::endl;

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

        // Compute the ordering permutation vector from the structural pattern of A
        m_solver.analyzePattern(m_A);
        // Compute the numerical factorization
        m_solver.factorize(m_A);
        if(m_solver.info()==Eigen::Success)
        {
            //Use the factors to solve the linear system
            qIter = m_solver.solve(m_b);
        }
        else
        {
            m_mesh.nodesList = nodesListSave;
            return false;
        }

        for(std::size_t n = 0 ; n < m_mesh.nodesList.size() ; ++n)
        {
            m_mesh.nodesList[n].states[0] = qIter[n];
            m_mesh.nodesList[n].states[1] = qIter[n + m_mesh.nodesList.size()];
            m_mesh.nodesList[n].states[2] = qIter[n + 2*m_mesh.nodesList.size()];
            m_mesh.nodesList[n].position[0] = nodesListSave[n].position[0]
                                            + m_currentDT*m_mesh.nodesList[n].states[0];
            m_mesh.nodesList[n].position[1] = nodesListSave[n].position[1]
                                            + m_currentDT*m_mesh.nodesList[n].states[1];
        }

        m_numIterSolverMsh++;

        if(m_numIterSolverMsh == 1)
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

        std::cout << " * Relative 2-norm of q: " << res << " vs "
                  << m_params.sideParams[0] << std::endl;

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

        if(m_numIterSolverMsh >= m_params.sideParams[1] || std::isnan(res))
        {
            m_mesh.nodesList = nodesListSave;
            return false;
        }
    }

    m_q = qIter;
    m_mesh.nodesList = nodesListSave;

    for(std::size_t n = 0 ; n < m_mesh.nodesList.size() ; ++n)
    {
        if(m_mesh.nodesList[n].isFree && !m_mesh.nodesList[n].isBound)
        {
            m_q(m_mesh.nodesList.size() + n) = m_q (m_mesh.nodesList.size() + n) -
                                                    m_currentDT*m_params.gravity*m_params.fluidParameters[0]*
                                                    m_params.hchar*m_params.hchar*0.5;
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
    std::vector<std::vector<Eigen::MatrixXd>> B;

    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        B.push_back(m_mesh.getB(elm));
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
    std::vector<Eigen::Triplet<double>> indexA;
    indexA.resize(99*m_mesh.nodesList.size());

    /********************************************************************************
                                    Build M/dt
    ********************************************************************************/
    m_M.resize(2*m_mesh.nodesList.size(), 2*m_mesh.nodesList.size());
    m_M.data().squeeze();
    index.reserve(2*3*3*m_mesh.getElementNumber());

    Eigen::MatrixXd Me(6, 6); Me.setZero();
    for (unsigned short k = 0 ; k < N.size() ; ++k)
        Me += N[k].transpose()*N[k]*GP2Dweight<double>[k];

    Me *= 0.5*(m_params.fluidParameters[0]/m_currentDT);

    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        for(unsigned short i = 0 ; i < 3 ; ++i)
        {
            for(unsigned short j = 0 ; j < 3 ; ++j)
            {
                index.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                       m_mesh.getElement(elm)[j],
                                                       Me(i,j)*m_mesh.getDetJ(elm)));

                index.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                       m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                       Me(i+3,j+3)*m_mesh.getDetJ(elm)));

                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                        m_mesh.getElement(elm)[j],
                                                        Me(i,j)*m_mesh.getDetJ(elm)));

                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                        m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                        Me(i+3,j+3)*m_mesh.getDetJ(elm)));
            }
        }
    }

    m_M.setFromTriplets(index.begin(), index.end());
    index.clear();

    /********************************************************************************
                                        Build K
    ********************************************************************************/
    m_K.resize(2*m_mesh.nodesList.size(), 2*m_mesh.nodesList.size());
    m_K.data().squeeze();
    index.reserve(6*6*m_mesh.getElementNumber());

    const Eigen::DiagonalMatrix<double, 3> ddev(2, 2, 1);

    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        Eigen::MatrixXd Ke = 0.5*m_params.fluidParameters[1]*B[elm][0].transpose()*ddev*B[elm][0]; //same matrices for all gauss point ^^

        for(unsigned short i = 0 ; i < 3 ; ++i)
        {
            for(unsigned short j = 0 ; j < 3 ; ++j)
            {
                index.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                       m_mesh.getElement(elm)[j],
                                                       Ke(i,j)*m_mesh.getDetJ(elm)));

                index.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                       m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                       Ke(i,j+3)*m_mesh.getDetJ(elm)));

                index.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                       m_mesh.getElement(elm)[j],
                                                       Ke(i+3,j)*m_mesh.getDetJ(elm)));

                index.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                       m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                       Ke(i+3,j+3)*m_mesh.getDetJ(elm)));

                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                        m_mesh.getElement(elm)[j],
                                                        Ke(i,j)*m_mesh.getDetJ(elm)));

                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                        m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                        Ke(i,j+3)*m_mesh.getDetJ(elm)));

                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                        m_mesh.getElement(elm)[j],
                                                        Ke(i+3,j)*m_mesh.getDetJ(elm)));

                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size(),
                                                        m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                        Ke(i+3,j+3)*m_mesh.getDetJ(elm)));
            }
        }
    }

    m_K.setFromTriplets(index.begin(), index.end());
    index.clear();

    /********************************************************************************
                                        Build D
    ********************************************************************************/
    m_D.resize(m_mesh.nodesList.size(), 2*m_mesh.nodesList.size());
    m_D.data().squeeze();
    index.reserve(3*6*m_mesh.getElementNumber());

    const Eigen::Vector3d m(1, 1, 0);

    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        Eigen::MatrixXd De(3,6); De.setZero();

        for (unsigned short k = 0 ; k < N.size() ; ++k)
            De += (N[k].topLeftCorner<1,3>()).transpose()*m.transpose()*B[elm][k]*GP2Dweight<double>[k];

        De *= 0.5;

        for(unsigned short i = 0 ; i < 3 ; ++i)
        {
            for(unsigned short j = 0 ; j < 3 ; ++j)
            {
                index.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                       m_mesh.getElement(elm)[j],
                                                       De(i,j)*m_mesh.getDetJ(elm)));

                index.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                       m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                       De(i,j+3)*m_mesh.getDetJ(elm)));

                //Part D of A
                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                        m_mesh.getElement(elm)[j],
                                                        De(i,j)*m_mesh.getDetJ(elm)));

                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                        m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                        De(i,j+3)*m_mesh.getDetJ(elm)));

                //Part -D^T of A
                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[j],
                                                        m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                        -De(i,j)*m_mesh.getDetJ(elm)));

                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                        m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                        -De(i,j+3)*m_mesh.getDetJ(elm)));
            }
        }
    }

    m_D.setFromTriplets(index.begin(), index.end());
    index.clear();


    /********************************************************************************
                                        Build C/dt
    ********************************************************************************/
    m_C.resize(m_mesh.nodesList.size(), 2*m_mesh.nodesList.size());
    m_C.data().squeeze();
    index.reserve(3*6*m_mesh.getElementNumber());

    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        Eigen::MatrixXd Ce(3,6); Ce.setZero();

        for (unsigned short k = 0 ; k < N.size() ; ++k)
            Ce +=(BleftP*B[elm][k]*BrightP).transpose()*N[k]*GP2Dweight<double>[k];

        Ce *= (0.5*m_tauPSPG[elm]/m_currentDT);

        for(unsigned short i = 0 ; i < 3 ; ++i)
        {
            for(unsigned short j = 0 ; j < 3 ; ++j)
            {
                index.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                       m_mesh.getElement(elm)[j],
                                                       Ce(i,j)*m_mesh.getDetJ(elm)));

                index.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                       m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                       Ce(i,j+3)*m_mesh.getDetJ(elm)));

                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                        m_mesh.getElement(elm)[j],
                                                        Ce(i,j)*m_mesh.getDetJ(elm)));

                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                        m_mesh.getElement(elm)[j] + m_mesh.nodesList.size(),
                                                        Ce(i,j+3)*m_mesh.getDetJ(elm)));
            }
        }
    }

    m_C.setFromTriplets(index.begin(), index.end());
    index.clear();

    /********************************************************************************
                                        Build L
    ********************************************************************************/
    m_L.resize(m_mesh.nodesList.size(), m_mesh.nodesList.size());
    m_L.data().squeeze();
    index.reserve(3*3*m_mesh.getElementNumber());

    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        Eigen::MatrixXd Le = 0.5*(m_tauPSPG[elm]/m_params.fluidParameters[0])
                             *(BleftP*B[elm][0]*BrightP).transpose()*(BleftP*B[elm][0]*BrightP);

        for(unsigned short i = 0 ; i < 3 ; ++i)
        {
            for(unsigned short j = 0 ; j < 3 ; ++j)
            {
                index.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i],
                                                       m_mesh.getElement(elm)[j],
                                                       Le(i,j)*m_mesh.getDetJ(elm)));

                indexA.push_back(Eigen::Triplet<double>(m_mesh.getElement(elm)[i] + 2*m_mesh.nodesList.size(),
                                                        m_mesh.getElement(elm)[j] + 2*m_mesh.nodesList.size(),
                                                        Le(i,j)*m_mesh.getDetJ(elm)));
            }
        }
    }

    m_L.setFromTriplets(index.begin(), index.end());
    index.clear();


    /********************************************************************************
                                        Build f
    ********************************************************************************/
    m_F.resize(2*m_mesh.nodesList.size());
    m_F.setZero();

    const Eigen::Vector2d b(0, -m_params.gravity);

    Eigen::MatrixXd Fe(6,1); Fe.setZero();

    for (unsigned short k = 0 ; k < N.size() ; ++k)
        Fe += N[k].transpose()*b*GP2Dweight<double>[k];

    Fe *= 0.5*m_params.fluidParameters[0];

    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        for(unsigned short i = 0 ; i < 3 ; ++i)
        {
            m_F(m_mesh.getElement(elm)[i]) += Fe(i)*m_mesh.getDetJ(elm);

            m_F(m_mesh.getElement(elm)[i] + m_mesh.nodesList.size()) +=
            Fe(i + 3)*m_mesh.getDetJ(elm);
        }
    }

    /********************************************************************************
                                        Build h
    ********************************************************************************/
    m_H.resize(m_mesh.nodesList.size());
    m_H.setZero();

    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        Eigen::VectorXd He = 0.5*m_tauPSPG[elm]*(BleftP*B[elm][0]*BrightP).transpose()*b;

        for(unsigned short i = 0 ; i < 3 ; ++i)
        {
            m_H(m_mesh.getElement(elm)[i]) += He(i)*m_mesh.getDetJ(elm);
        }
    }

    /********************************************************************************
                                        Compute A and b
    ********************************************************************************/
    m_A.resize(3*m_mesh.nodesList.size(), 3*m_mesh.nodesList.size());
    m_A.data().squeeze();
    m_A.setFromTriplets(indexA.begin(), indexA.end());

    m_b.resize(3*m_mesh.nodesList.size());
    m_b.setZero();
    m_b << m_F + m_M*m_qprev.head(2*m_mesh.nodesList.size()), m_H + m_C*m_qprev.head(2*m_mesh.nodesList.size());
}

void Solver::_computeTauPSPG()
{
    m_tauPSPG.resize(m_mesh.getElementNumber());

    double rho = m_params.fluidParameters[0];
    double mu = m_params.fluidParameters[1];

    double U = 0;
    std::size_t trueNnodes = 0;
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

    for(std::size_t elm = 0 ; elm < m_mesh.getElementNumber() ; ++elm)
    {
        double h = std::sqrt(2*m_mesh.getDetJ(elm)/M_PI);

        m_tauPSPG[elm] = 1/std::sqrt((2/m_currentDT)*(2/m_currentDT)
                                 + (2*U/h)*(2*U/h)
                                 + 9*(4*mu/(h*h*rho))*(4*mu/(h*h*rho)));
    }
}

std::vector<bool> Solver::_getIndices() const
{
    std::vector<bool> indices(3*m_mesh.nodesList.size());

    for (std::size_t i = 0 ; i < m_mesh.nodesList.size() ; ++i)
    {
        if(m_mesh.nodesList[i].isBound || m_mesh.nodesList[i].isFree)
        {
            indices[i] = true;
            indices[i + m_mesh.nodesList.size()] = true;
        }
        else
        {
            indices[i] = false;
            indices[i + m_mesh.nodesList.size()] = false;
        }

        if(m_mesh.nodesList[i].isFree)
            indices[i + 2*m_mesh.nodesList.size()] = true;
        else
            indices[i + 2*m_mesh.nodesList.size()] = false;
    }

    return indices;
}
