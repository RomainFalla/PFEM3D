#include <chrono>
#include <iomanip>
#include <iostream>

#include <gmsh.h>

#include "SolverCompressible.hpp"

using Clock = std::chrono::high_resolution_clock;
using TimeType = std::chrono::time_point<std::chrono::high_resolution_clock>;

static void displayDT(TimeType startTime, TimeType endTime, std::string text)
{
    auto ellapsedTimeMeasure = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << text << static_cast<double>(ellapsedTimeMeasure.count())/1000.0 << " s" << std::endl;
}


SolverCompressible::SolverCompressible(const SolverCompCreateInfo& solverCompInfos) :
Solver(solverCompInfos.solverInfos),
m_rho0(solverCompInfos.rho0),
m_mu(solverCompInfos.mu),
m_K0(solverCompInfos.K0),
m_K0prime(solverCompInfos.K0prime),
m_pInfty(solverCompInfos.pInfty),
m_securityCoeff(solverCompInfos.securityCoeff),
m_strongContinuity(solverCompInfos.strongContinuity),
m_nextTimeToRemesh(solverCompInfos.solverInfos.maxDT)
{
    m_solverType = SOLVER_TYPE::WeaklyCompressible;
    unsigned short dim = m_mesh.getDim();

    m_statesNumber = 2*dim + 2;

    m_mesh.setStatesNumber(m_statesNumber);

    m_lua["rho0"] = m_rho0;
    m_lua["K0"] = m_K0;
    m_lua["K0prime"] = m_K0prime;
    m_lua["mu"] = m_mu;
    m_lua["pInfty"] = m_pInfty;

    m_MrhoPrev.resize(dim + 1, dim + 1); m_MrhoPrev.setZero();
    for(unsigned short k = 0 ; k < m_N2.size() ; ++k)
    {
        m_MrhoPrev += (m_N2[k].topLeftCorner(1, dim + 1)).transpose()*m_N2[k].topLeftCorner(1,  dim + 1)*m_w2[k];
    }
    m_MrhoPrev *= m_mesh.getRefElementSize(dim);

    m_MrhoLumpedPrev.resize(dim + 1); m_MrhoLumpedPrev.setZero();
    for(unsigned short i = 0 ; i < m_MrhoPrev.rows() ; ++i)
    {
        for(unsigned short k = 0 ; k < m_MrhoPrev.cols() ; ++k)
        {
                m_MrhoLumpedPrev.diagonal()[i] += m_MrhoPrev(i, k);
        }
    }

    if(dim == 2)
    {
        m_ddev.resize(3, 3);
        m_ddev <<  4.0/3, -2.0/3, 0,
                  -2.0/3,  4.0/3, 0,
                       0,      0, 1;
    }
    else
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

    setInitialCondition();

    m_mesh.setComputeNormalCurvature(false);
}

void SolverCompressible::applyBoundaryConditionsCont(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invMrho, Eigen::VectorXd& Frho)
{
    assert(m_mesh.getNodesCount() != 0);

    auto& invMrhoDiag = invMrho.diagonal();

    //Do not parallelize this
    for (std::size_t n = 0 ; n < m_mesh.getNodesCount() ; ++n)
    {
        const Node& node = m_mesh.getNode(n);

        if(node.isFree())
        {
            Frho(n) = m_rho0;

            invMrhoDiag[n] = 1;
        }
        else if(m_strongPAtFS && node.isOnFreeSurface())
        {
            Frho(n) = m_rho0;

            invMrhoDiag[n] = 1;
        }
    }
}

void SolverCompressible::applyBoundaryConditionsMom(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F)
{
    assert(m_mesh.getNodesCount() != 0);

    const uint8_t dim = m_mesh.getDim();
    const std::size_t nodesCount = m_mesh.getNodesCount();

    auto& invMDiag = invM.diagonal();

    //Do not parallelize this
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_mesh.getNode(n);

        if(node.isFree() && !node.isBound())
        {
            for(unsigned short d = 0 ; d < dim - 1 ; ++d)
            {
                F(n + d*nodesCount) = 0;
                invMDiag[n + d*nodesCount] = 1;
            }

             F(n + (dim - 1)*nodesCount) = - m_gravity;
             invMDiag[n + (dim - 1)*nodesCount] = 1;
        }
        else if(node.isBound())
        {
            std::vector<double> result;
            result = m_lua[m_mesh.getNodeType(n)](node.getPosition(),
                                                  m_mesh.getBoundNodeInitPos(n),
                                                  m_currentTime + m_currentDT).get<std::vector<double>>();

            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                F(n + d*nodesCount) = result[d];
                invMDiag[n + d*nodesCount] = 1;
            }
        }
    }
}

void SolverCompressible::displaySolverParams() const noexcept
{
    std::cout << "Initial nodes number: " << m_mesh.getNodesCount() << "\n";
    std::cout << "Initial elements number: " << m_mesh.getElementsCount() << "\n";
    m_mesh.displayToConsole();
#if defined(_OPENMP)
    std::cout << "Number of OpenMP threads: " << m_numOMPThreads << "/" << omp_get_num_procs() << "\n";
#endif
    std::cout << "Gravity: " << m_gravity << " m/s^2" << "\n";
    std::cout << "Density: " << m_rho0 << " kg/m^3" << "\n";
    std::cout << "Viscosity: " << m_mu << " Pa s" << "\n";
    std::cout << "K0: " << m_K0 << " Pa" << "\n";
    std::cout << "K0prime: " << m_K0prime << "\n";
    std::cout << "pInfty: " << m_pInfty << "\n";
    std::cout << "End simulation time: " << m_endTime << " s" << "\n";
    std::cout << "Adapt time step: " << (m_adaptDT ? "yes" : "no") << "\n";
    if(m_adaptDT)
    {
        std::cout << "Maximum time step: " << m_maxDT << " s" << "\n";
        std::cout << "Security coefficient: " << m_securityCoeff << "\n";
    }
    else
        std::cout << "Time step: " << m_maxDT << " s" << "\n";

    std::cout << std::flush;
}

void SolverCompressible::solveProblem(bool verboseOutput)
{
    std::cout   << "================================================================"
                << "\n"
                << "               EXECUTING PFEM COMPRESSIBLE SOLVER               "
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
            //if(m_currentStep%100 == 0 || m_currentTime + m_currentDT >= m_endTime)
            //{
                std::cout << std::fixed << std::setprecision(3);
                std::cout << "\r" << "Solving time step: " << m_currentTime + m_currentDT
                          << "/" << m_endTime << " s, dt = ";
                std::cout << std::scientific;
                std::cout << m_currentDT << " s" << std::flush;
            //}

        }

        solveCurrentTimeStep(verboseOutput);

        for(auto& extractor : m_pExtractor)
        {
            extractor->update();
        }

        if(m_adaptDT)
        {
            double U = 0;
            double c = 0;
            for(std::size_t n = 0 ; n < m_mesh.getNodesCount() ; ++n)
            {
                const Node& node = m_mesh.getNode(n);
                if(node.isBound() && node.isFree())
                    continue;

                double localU = 0;
                for(uint16_t d = 0 ; d < m_mesh.getDim() ; ++d)
                    localU += node.getState(d)*node.getState(d);

                U = std::max(localU, U);

                c = std::max((m_K0 + m_K0prime*node.getState(m_mesh.getDim()))/node.getState(m_mesh.getDim() + 1), c);
            }
            U = std::sqrt(U);
            c = std::sqrt(c);

            m_currentDT = std::min(m_maxDT, m_securityCoeff*m_mesh.getHchar()/std::max(U, c));
        }
    }

    std::cout << std::endl;
}

bool SolverCompressible::solveCurrentTimeStep(bool verboseOutput)
{
    TimeType startTime, endTime, startTimeMeasure, endTimeMeasure;
    startTime = Clock::now();

    const unsigned short dim = m_mesh.getDim();

    Eigen::VectorXd qVPrev = getQFromNodesStates(0, dim - 1);           //The precedent speed.
    Eigen::VectorXd qAccPrev = getQFromNodesStates(dim + 2, 2*dim + 1);  //The precedent acceleration.

    Eigen::VectorXd Frho;                                 //The rhos of the continuity equation.

    if(m_strongContinuity)
        buildFrho(Frho);

    Eigen::VectorXd qV1half = qVPrev + 0.5*m_currentDT*qAccPrev;

    {
        setNodesStatesfromQ(qV1half, 0, dim - 1);
        Eigen::VectorXd deltaPos = qV1half*m_currentDT;
        m_mesh.updateNodesPosition(std::vector<double> (deltaPos.data(), deltaPos.data() + deltaPos.cols()*deltaPos.rows()));
    }

    {
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> invMrho; //The mass matrix of the continuity.
        buildMatricesCont(invMrho, Frho);
        applyBoundaryConditionsCont(invMrho, Frho);

        Eigen::VectorXd qRho(m_mesh.getNodesCount());
        qRho = invMrho*Frho;
        setNodesStatesfromQ(qRho, dim + 1, dim + 1);
        Eigen::VectorXd qP = getPFromRhoTaitMurnagham(qRho);

        setNodesStatesfromQ(qP, dim, dim);
    }

    {
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> invM; //The mass matrix for momentum equation.
        Eigen::VectorXd F;                                  //The rhs of the momentum equation.

        buildMatricesMom(invM, F);
        applyBoundaryConditionsMom(invM, F);

        qAccPrev = invM*F;

        qVPrev = qV1half + 0.5*m_currentDT*qAccPrev;

        setNodesStatesfromQ(qVPrev, 0, dim - 1);
        setNodesStatesfromQ(qAccPrev, dim + 2, 2*dim + 1);
    }

    m_currentTime += m_currentDT;
    m_currentStep++;

    endTime = Clock::now();
    if(verboseOutput)
        displayDT(startTime, endTime, "Problem solved in ");

    if(m_currentTime > m_nextTimeToRemesh)
    {
        startTime = Clock::now();
        m_mesh.remesh(verboseOutput);
        endTime = Clock::now();
        if(verboseOutput)
            displayDT(startTime, endTime, "Remeshing done in ");

        m_nextTimeToRemesh += m_maxDT;
    }

    return true;
}

void SolverCompressible::buildFrho(Eigen::VectorXd& Frho)
{
    const unsigned short dim = m_mesh.getDim();

    Frho.resize(m_mesh.getNodesCount()); Frho.setZero();

    std::vector<Eigen::VectorXd> Frhoe(m_mesh.getElementsCount());

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_mesh.getElementsCount() ; ++elm)
    {
        const Element& element = m_mesh.getElement(elm);

        Eigen::VectorXd Rho = getElementState(element, dim + 1);

        Eigen::MatrixXd Mrhoe = m_MrhoPrev*element.getDetJ();

        Frhoe[elm] = Mrhoe*Rho;
    }
    Eigen::setNbThreads(m_numOMPThreads);

    for(std::size_t elm = 0 ; elm < m_mesh.getElementsCount() ; ++elm)
    {
        const Element& element = m_mesh.getElement(elm);

        for(uint8_t i = 0 ; i < dim + 1 ; ++i)
            Frho(element.getNodeIndex(i)) += Frhoe[elm](i);
    }
}

void SolverCompressible::buildMatricesCont(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invMrho, Eigen::VectorXd& Frho)
{
    const unsigned short dim = m_mesh.getDim();

    invMrho.resize(m_mesh.getNodesCount()); invMrho.setZero();

    std::vector<Eigen::DiagonalMatrix<double,Eigen::Dynamic>> MrhoeLumped(m_mesh.getElementsCount());
    std::vector<Eigen::VectorXd> Frhoe;

    if(!m_strongContinuity)
    {
        Frho.resize(m_mesh.getNodesCount()); Frho.setZero();
        Frhoe.resize((m_mesh.getElementsCount()));
    }

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_mesh.getElementsCount() ; ++elm)
    {
        const Element& element = m_mesh.getElement(elm);

        MrhoeLumped[elm] = element.getDetJ()*m_MrhoLumpedPrev;

        if(!m_strongContinuity)
        {
            Eigen::VectorXd Rho = getElementState(element, dim + 1);

            Eigen::VectorXd V((dim + 1)*dim);
            if(dim == 2)
            {
                V << getElementState(element, 0),
                     getElementState(element, 1);
            }
            else
            {
                 V << getElementState(element, 0),
                      getElementState(element, 1),
                      getElementState(element, 2);
            }

            Eigen::MatrixXd Be = getB(element);

            //De = S rho Np^T m Bv dV
            Eigen::MatrixXd Drhoe(dim + 1, dim*(dim + 1)); Drhoe.setZero();

            for (unsigned short k = 0 ; k < m_N2.size() ; ++k)
            {
                double rho = (m_N2[k].topLeftCorner(1, dim + 1)*Rho).value();
                Drhoe += rho*m_N2[k].topLeftCorner(1, dim + 1).transpose()*m_m.transpose()*Be*m_w2[k];
            }

            Drhoe *= m_mesh.getRefElementSize(dim)*element.getDetJ();

            Frhoe[elm] = element.getDetJ()*m_MrhoPrev*Rho - m_currentDT*Drhoe*V;
        }
    }
    Eigen::setNbThreads(m_numOMPThreads);

    auto& invMrhoDiag = invMrho.diagonal();

    for(std::size_t elm = 0 ; elm < m_mesh.getElementsCount() ; ++elm)
    {
        const Element& element = m_mesh.getElement(elm);

        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
        {
            invMrhoDiag[element.getNodeIndex(i)] += MrhoeLumped[elm].diagonal()[i];

            if(!m_strongContinuity)
            {
                Frho(element.getNodeIndex(i)) += Frhoe[elm](i);
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(auto i = 0 ; i < invMrho.rows() ; ++i)
         invMrhoDiag[i] = 1/invMrhoDiag[i];
}

void SolverCompressible::buildMatricesMom(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F)
{
    const unsigned short dim = m_mesh.getDim();
    const std::size_t elementsCount = m_mesh.getElementsCount();
    const std::size_t nodesCount = m_mesh.getNodesCount();

    invM.resize(dim*nodesCount); invM.setZero();
    std::vector<Eigen::MatrixXd> Me(elementsCount);

    F.resize(dim*nodesCount); F.setZero();
    std::vector<Eigen::VectorXd> FTote(elementsCount);

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < elementsCount ; ++elm)
    {
        const Element& element = m_mesh.getElement(elm);

        Eigen::VectorXd Rho = getElementState(element, dim + 1);

        Eigen::VectorXd V((dim + 1)*dim);
        if(dim == 2)
        {
            V << getElementState(element, 0),
                 getElementState(element, 1);
        }
        else
        {
             V << getElementState(element, 0),
                  getElementState(element, 1),
                  getElementState(element, 2);
        }

        Eigen::VectorXd P = getElementState(element, dim);

        Eigen::MatrixXd Be = getB(element);

        Me[elm].resize(dim*(dim + 1), dim*(dim + 1)); Me[elm].setZero();

        for(unsigned short k = 0 ; k < m_N3.size() ; ++k)
        {
            double rho = (m_N3[k].topLeftCorner(1, dim + 1)*Rho).value();
            Me[elm] += rho*m_N3[k].transpose()*m_N3[k]*m_w3[k];
        }

        Me[elm] *= m_mesh.getRefElementSize(dim)*element.getDetJ();

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
        Eigen::MatrixXd Ke = m_mesh.getRefElementSize(dim)*Be.transpose()*m_ddev*Be*element.getDetJ(); //same matrices for all gauss point ^^

        //De = S Np^T m Bv dV
        Eigen::MatrixXd De(dim + 1, dim*(dim + 1)); De.setZero();

        for (unsigned short k = 0 ; k < m_N1.size() ; ++k)
            De += (m_N1[k].topLeftCorner(1, dim + 1)).transpose()*m_m.transpose()*Be*m_w1[k];

        De *= m_mesh.getRefElementSize(dim)*element.getDetJ();

        //Fe = S Nv^T bodyforce dV
        Eigen::MatrixXd Fe(dim*(dim + 1), 1); Fe.setZero();
        //Currently body force is constant => Fe is of order 2. Might need to chage this if one use a user-def body force
        for (unsigned short k = 0 ; k < m_N2.size() ; ++k)
        {
            double rho = (m_N2[k].topLeftCorner(1, dim + 1)*Rho).value();
            Fe += rho*m_N2[k].transpose()*m_bodyForces*m_w2[k];
        }

        Fe *= m_mesh.getRefElementSize(dim)*element.getDetJ();

        FTote[elm] = -Ke*V + De.transpose()*P + Fe;
    }
    Eigen::setNbThreads(m_numOMPThreads);

    auto& invMDiag = invM.diagonal();

    for(std::size_t elm = 0 ; elm < elementsCount ; ++elm)
    {
        const Element& element = m_mesh.getElement(elm);

        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                /********************************************************************
                                             Build M
                ********************************************************************/
                invMDiag[element.getNodeIndex(i) + d*nodesCount] += Me[elm](i + d*(dim + 1), i + d*(dim + 1));

                /************************************************************************
                                                Build f
                ************************************************************************/
                F(element.getNodeIndex(i) + d*nodesCount) += FTote[elm](i + d*(dim + 1));
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(std::size_t i = 0 ; i < dim*nodesCount ; ++i)
        invMDiag[i] = 1/invMDiag[i];
}
