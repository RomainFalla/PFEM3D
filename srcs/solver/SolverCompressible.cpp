#include <iomanip>
#include <iostream>

#include <gmsh.h>

#include "SolverCompressible.hpp"


SolverCompressible::SolverCompressible(const nlohmann::json& j, const std::string& mshName, const std::string& resultsName) :
Solver(j, mshName, resultsName)
{
    unsigned short dim = m_mesh.getMeshDim();

    m_statesNumber = 2*dim + 2;

    m_mesh.setStatesNumber(m_statesNumber);

    m_strongContinuity        = j["Solver"]["strongContinuity"].get<bool>();

    m_rho0              = j["Solver"]["Fluid"]["rho0"].get<double>();
    m_mu                = j["Solver"]["Fluid"]["mu"].get<double>();
    m_K0                = j["Solver"]["Fluid"]["K0"].get<double>();
    m_K0prime           = j["Solver"]["Fluid"]["K0prime"].get<double>();
    m_pInfty            = j["Solver"]["Fluid"]["pInfty"].get<double>();

    m_securityCoeff      = j["Solver"]["Time"]["securityCoeff"].get<double>();

    if(dim == 2)
        m_whatCanBeWriten = {"u", "v", "p", "rho", "ax", "ay", "ke", "velocity"};
    else if(dim == 3)
        m_whatCanBeWriten = {"u", "v", "w", "p", "rho", "ax", "ay", "az", "ke", "velocity"};

    std::vector<std::string> whatToWrite = j["Solver"]["whatToWrite"];
    m_whatToWrite.resize(2*dim + 4, false);

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

    m_sumNTN.resize(dim + 1, dim + 1); m_sumNTN.setZero();
    for(unsigned short k = 0 ; k < m_N.size() ; ++k)
    {
        for(unsigned short k = 0 ; k < m_N.size() ; ++k)
            m_sumNTN += (m_N[k].topLeftCorner(1, dim + 1)).transpose()*m_N[k].topLeftCorner(1,  dim + 1)*m_mesh.getGaussWeight(k);
    }

    if(dim == 2)
    {
        m_ddev.resize(3, 3);
        m_ddev <<  4.0/3, -2.0/3, 0,
                  -2.0/3,  4.0/3, 0,
                       0,      0, 1;
    }
    else if(dim == 3)
    {
        m_ddev.resize(6, 6);
        m_ddev <<  4.0/3, -2.0/3, -2.0/3, 0, 0, 0,
                  -2.0/3,  4.0/3, -2.0/3, 0, 0, 0,
                  -2.0/3, -2.0/3,  4.0/3, 0, 0, 0,
                       0,      0,      0, 1, 0, 0,
                       0,      0,      0, 0, 1, 0,
                       0,      0,      0, 0, 0, 1;
    }

    m_ddev *= m_mu;
}

SolverCompressible::~SolverCompressible()
{
    gmsh::finalize();
}

void SolverCompressible::applyBoundaryConditionsMom()
{
    assert(m_mesh.getNodesNumber() != 0);

    const unsigned short dim = m_mesh.getMeshDim();

    //Do not parallelize this
    for (std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        if(m_mesh.isNodeFree(n) && !m_mesh.isNodeBound(n))
        {
            for(unsigned short d = 0 ; d < dim - 1 ; ++d)
            {
                m_F(n + d*m_mesh.getNodesNumber()) = 0;
                m_invM.diagonal()[n + d*m_mesh.getNodesNumber()] = 1;
            }

             m_F(n + (dim - 1)*m_mesh.getNodesNumber()) = - m_gravity;
             m_invM.diagonal()[n + (dim - 1)*m_mesh.getNodesNumber()] = 1;
        }
        else if(m_mesh.isNodeBound(n))
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                m_F(n + d*m_mesh.getNodesNumber()) = 0;
                m_invM.diagonal()[n + d*m_mesh.getNodesNumber()] = 1;
            }
        }
    }
}

void SolverCompressible::applyBoundaryConditionsCont()
{
    assert(m_mesh.getNodesNumber() != 0);

    //Do not parallelize this
    for (std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        if(m_mesh.isNodeFree(n))
        {
           m_Frho(n) = m_mesh.getNodeState(n, m_mesh.getMeshDim() + 1);

           m_invMrho.diagonal()[n] = 1;
        }
    }
}

void SolverCompressible::displaySolverParams() const
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
    std::cout << "Density: " << m_rho0 << " kg/m^3" << std::endl;
    std::cout << "Viscosity: " << m_mu << " Pa s" << std::endl;
    std::cout << "K0: " << m_K0 << " Pa" << std::endl;
    std::cout << "K0prime: " << m_K0prime << std::endl;
    std::cout << "pInfty: " << m_pInfty << std::endl;
    std::cout << "End simulation time: " << m_endTime << " s" << std::endl;
    std::cout << "Write solution every: " << m_timeBetweenWriting << " s" << std::endl;
    std::cout << "Adapt time step: " << (m_adaptDT ? "yes" : "no") << std::endl;
    if(m_adaptDT)
    {
        std::cout << "Maximum time step: " << m_maxDT << " s" << std::endl;
        std::cout << "Security coefficient: " << m_securityCoeff << std::endl;
    }
    else
        std::cout << "Time step: " << m_maxDT << " s" << std::endl;
}

void SolverCompressible::setInitialCondition()
{
    assert(m_mesh.getNodesNumber() != 0);

    const unsigned short dim = m_mesh.getMeshDim();

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        for(unsigned short d = 0 ; d < dim + 1 ; ++d)
        {
            if(!m_mesh.isNodeBound(n) || m_mesh.isNodeFluidInput(n))
                m_mesh.setNodeState(n, d, m_initialCondition[d]);
            else
                m_mesh.setNodeState(n, d, 0);
        }

        m_mesh.setNodeState(n, dim + 1, m_initialCondition[dim + 1]);
        for(unsigned short d = 0 ; d < dim ; ++d)
        {
            m_mesh.setNodeState(n, dim + 2 + d, 0);
        }
    }

    m_qVPrev = getQFromNodesStates(0, dim - 1);
    m_qAccPrev = getQFromNodesStates(dim + 2, 2*dim + 1);
}

void SolverCompressible::solveProblem()
{
    std::cout   << "================================================================"
                << std::endl
                << "               EXECUTING PFEM COMPRESSIBLE SOLVER               "
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

        solveCurrentTimeStep();

        if(m_currentTime >= m_nextWriteTrigger)
        {
            writeData();
            m_nextWriteTrigger += m_timeBetweenWriting;
        }

        if(m_adaptDT)
        {
            double U = 0;
            double c = 0;
            for(std::size_t n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
            {
                if(m_mesh.isNodeBound(n) && m_mesh.isNodeFree(n))
                    continue;

                double localU = 0;
                for(unsigned short d = 0 ; d < m_mesh.getMeshDim() ; ++d)
                    localU += m_mesh.getNodeState(n, d)*m_mesh.getNodeState(n, d);

                U = std::max(localU, U);

                c = std::max((m_K0 + m_K0prime *m_mesh.getNodeState(n, m_mesh.getMeshDim()))/m_mesh.getNodeState(n, m_mesh.getMeshDim() + 1), c);
            }
            U = std::sqrt(U);
            c = std::sqrt(c);

            m_currentDT = std::min(m_maxDT, m_securityCoeff*m_mesh.getHchar()/std::max(U, c));
        }
    }

    std::cout << std::endl;
}

bool SolverCompressible::solveCurrentTimeStep()
{
    assert(m_qVPrev.size() == m_qAccPrev.size());
    assert(m_qVPrev.size() == m_mesh.getMeshDim()*m_mesh.getNodesNumber());

    const unsigned short dim = m_mesh.getMeshDim();

    if(m_strongContinuity)
        buildFrho();

    Eigen::VectorXd qV1half = m_qVPrev + 0.5*m_currentDT*m_qAccPrev;

    setNodesStatesfromQ(qV1half, 0, dim - 1);
    Eigen::VectorXd deltaPos = qV1half*m_currentDT;
    m_mesh.updateNodesPosition(std::vector<double> (deltaPos.data(), deltaPos.data() + deltaPos.cols()*deltaPos.rows()));

    buildMatricesCont();
    applyBoundaryConditionsCont();

    Eigen::VectorXd qRho(m_mesh.getNodesNumber());
    qRho = m_invMrho*m_Frho;
    setNodesStatesfromQ(qRho, dim + 1, dim + 1);
    Eigen::VectorXd qP = getPFromRhoTaitMurnagham(qRho);

    setNodesStatesfromQ(qP, dim, dim);

    buildMatricesMom();
    applyBoundaryConditionsMom();

    m_qAccPrev = m_invM*m_F;

    m_qVPrev = qV1half + 0.5*m_currentDT*m_qAccPrev;

    setNodesStatesfromQ(m_qVPrev, 0, dim - 1);
    setNodesStatesfromQ(m_qAccPrev, dim + 2, 2*dim + 1);

    m_currentTime += m_currentDT;
    m_currentStep++;

    //Remeshing step
    m_mesh.remesh();

    //We have to compute qPrev here due to new nodes !
    m_qVPrev = getQFromNodesStates(0, dim - 1);
    m_qAccPrev = getQFromNodesStates(dim + 2, 2*dim + 1);

    return true;
}

void SolverCompressible::buildFrho()
{
    const unsigned short dim = m_mesh.getMeshDim();

    m_Frho.resize(m_mesh.getNodesNumber()); m_Frho.setZero();

    std::vector<Eigen::VectorXd> Frhoe(m_mesh.getElementsNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        Eigen::VectorXd Rho = getElementState(elm, dim + 1);

        Eigen::MatrixXd Mrhoe = m_mesh.getRefElementSize()*m_mesh.getElementDetJ(elm)*m_sumNTN;

        Frhoe[elm] = Mrhoe*Rho;
    }

    for(std::size_t elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
            m_Frho(m_mesh.getElement(elm)[i]) += Frhoe[elm](i);
    }
}

void SolverCompressible::buildMatricesCont()
{
    const unsigned short dim = m_mesh.getMeshDim();

    m_invMrho.resize(m_mesh.getNodesNumber()); m_invMrho.setZero();

    std::vector<Eigen::MatrixXd> Mrhoe(m_mesh.getElementsNumber());
    std::vector<Eigen::VectorXd> Frhoe;

    if(!m_strongContinuity)
    {
        m_Frho.resize(m_mesh.getNodesNumber()); m_Frho.setZero();
        Frhoe.resize((m_mesh.getElementsNumber()));
    }

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        Mrhoe[elm] = m_mesh.getRefElementSize()*m_mesh.getElementDetJ(elm)*m_sumNTN;

        if(!m_strongContinuity)
        {
            Eigen::VectorXd Rho = getElementState(elm, dim + 1);

            Eigen::VectorXd V((dim + 1)*dim);
            if(dim == 2)
            {
                V << getElementState(elm, 0),
                     getElementState(elm, 1);
            }
            else if(dim == 3)
            {
                 V << getElementState(elm, 0),
                      getElementState(elm, 1),
                      getElementState(elm, 2);
            }

            Eigen::MatrixXd Be = getB(elm);

            //De = S rho Np^T m Bv dV
            Eigen::MatrixXd Drhoe(dim + 1, dim*(dim + 1)); Drhoe.setZero();

            for (unsigned short k = 0 ; k < m_N.size() ; ++k)
            {
                double rho = (m_N[k].topLeftCorner(1, dim + 1)*Rho).value();
                Drhoe += rho*m_N[k].topLeftCorner(1, dim + 1).transpose()*m_m.transpose()*Be*m_mesh.getGaussWeight(k);
            }

            Drhoe *= m_mesh.getRefElementSize()*m_mesh.getElementDetJ(elm);

            Frhoe[elm] = Mrhoe[elm]*Rho - m_currentDT*Drhoe*V;
        }

        for(unsigned short i = 0 ; i < Mrhoe[elm].rows() ; ++i)
        {
            for(unsigned short j = 0 ; j < Mrhoe[elm].cols() ; ++j)
            {
                if(i != j)
                {
                    Mrhoe[elm](i, i) += Mrhoe[elm](i, j);
                    Mrhoe[elm](i, j) = 0;
                }
            }
        }
    }

    for(std::size_t elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
        {
            m_invMrho.diagonal()[m_mesh.getElement(elm)[i]] += Mrhoe[elm](i,i);

            if(!m_strongContinuity)
            {
                m_Frho(m_mesh.getElement(elm)[i]) += Frhoe[elm](i);
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(std::size_t i = 0 ; i < m_invMrho.rows() ; ++i)
         m_invMrho.diagonal()[i] = 1/m_invMrho.diagonal()[i];
}

void SolverCompressible::buildMatricesMom()
{
    const unsigned short dim = m_mesh.getMeshDim();

    m_invM.resize(dim*m_mesh.getNodesNumber()); m_invM.setZero();
    std::vector<Eigen::MatrixXd> Me(m_mesh.getElementsNumber());

    m_F.resize(dim*m_mesh.getNodesNumber()); m_F.setZero();
    std::vector<Eigen::VectorXd> FTote(m_mesh.getElementsNumber());

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        Eigen::VectorXd Rho = getElementState(elm, dim + 1);

        Eigen::VectorXd V((dim + 1)*dim);
        if(dim == 2)
        {
            V << getElementState(elm, 0),
                 getElementState(elm, 1);
        }
        else if(dim == 3)
        {
             V << getElementState(elm, 0),
                  getElementState(elm, 1),
                  getElementState(elm, 2);
        }

        Eigen::VectorXd P = getElementState(elm, dim);

        Eigen::MatrixXd Be = getB(elm);

        Me[elm].resize(dim*(dim + 1), dim*(dim + 1)); Me[elm].setZero();

        for(unsigned short k = 0 ; k < m_N.size() ; ++k)
        {
            double rho = (m_N[k].topLeftCorner(1, dim + 1)*Rho).value();
            Me[elm] += rho*m_N[k].transpose()*m_N[k]*m_mesh.getGaussWeight(k);
        }

        Me[elm] *= m_mesh.getRefElementSize()*m_mesh.getElementDetJ(elm);

        for(unsigned short i = 0 ; i < Me[elm].rows() ; ++i)
        {
            for(unsigned short j = 0 ; j < Me[elm].cols() ; ++j)
            {
                if(i != j)
                {
                    Me[elm](i, i) += Me[elm](i, j);
                    Me[elm](i, j) = 0;
                }
            }
        }

        //Ke = S Bv^T ddev Bv dV
        Eigen::MatrixXd Ke = m_mesh.getRefElementSize()*Be.transpose()*m_ddev*Be*m_mesh.getElementDetJ(elm); //same matrices for all gauss point ^^

        //De = S Np^T m Bv dV
        Eigen::MatrixXd De(dim + 1, dim*(dim + 1)); De.setZero();

        for (unsigned short k = 0 ; k < m_N.size() ; ++k)
            De += (m_N[k].topLeftCorner(1, dim + 1)).transpose()*m_m.transpose()*Be*m_mesh.getGaussWeight(k);

        De *= m_mesh.getRefElementSize()*m_mesh.getElementDetJ(elm);

        //Fe = S Nv^T bodyforce dV
        Eigen::MatrixXd Fe(dim*(dim + 1), 1); Fe.setZero();

        for (unsigned short k = 0 ; k < m_N.size() ; ++k)
        {
            double rho = (m_N[k].topLeftCorner(1, dim + 1)*Rho).value();
            Fe += rho*m_N[k].transpose()*m_bodyForces*m_mesh.getGaussWeight(k);
        }

        Fe *= m_mesh.getRefElementSize()*m_mesh.getElementDetJ(elm);

        Eigen::MatrixXd FeTot(dim*(dim + 1), 1); FeTot.setZero();
        FTote[elm] = -Ke*V + De.transpose()*P + Fe;
    }

    for(std::size_t elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                /********************************************************************
                                             Build M
                ********************************************************************/
                m_invM.diagonal()[m_mesh.getElement(elm)[i] + d*m_mesh.getNodesNumber()] += Me[elm](i + d*(dim + 1), i + d*(dim + 1));

                /************************************************************************
                                                Build f
                ************************************************************************/
                m_F(m_mesh.getElement(elm)[i] + d*m_mesh.getNodesNumber()) += FTote[elm](i + d*(dim + 1));

            }
        }
    }

    #pragma omp parallel for default(shared)
    for(std::size_t i = 0 ; i < dim*m_mesh.getNodesNumber() ; ++i)
        m_invM.diagonal()[i] = 1/m_invM.diagonal()[i];
}
