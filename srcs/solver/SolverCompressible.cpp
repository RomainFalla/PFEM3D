#include <chrono>
#include <iomanip>
#include <iostream>

#include <gmsh.h>

#include "SolverCompressible.hpp"


SolverCompressible::SolverCompressible(const nlohmann::json& j, const std::string& mshName) :
Solver(j, mshName)
{
    m_solverType = WeaklyCompressible;
    unsigned short dim = m_mesh.getDim();

    m_statesNumber = 2*dim + 2;

    m_mesh.setStatesNumber(m_statesNumber);

    m_strongContinuity        = j["Solver"]["strongContinuity"].get<bool>();

    m_rho0              = j["Solver"]["Fluid"]["rho0"].get<double>();
    m_mu                = j["Solver"]["Fluid"]["mu"].get<double>();
    m_K0                = j["Solver"]["Fluid"]["K0"].get<double>();
    m_K0prime           = j["Solver"]["Fluid"]["K0prime"].get<double>();
    m_pInfty            = j["Solver"]["Fluid"]["pInfty"].get<double>();

    m_securityCoeff      = j["Solver"]["Time"]["securityCoeff"].get<double>();

    std::vector<std::string> whatCanBeWritten;
    if(dim == 2)
        whatCanBeWritten = {"u", "v", "p", "rho", "ax", "ay", "ke", "velocity"};
    else
        whatCanBeWritten = {"u", "v", "w", "p", "rho", "ax", "ay", "az", "ke", "velocity"};

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
}

SolverCompressible::~SolverCompressible()
{
}

void SolverCompressible::applyBoundaryConditionsCont(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invMrho, Eigen::VectorXd& Frho)
{
    assert(m_mesh.getNodesNumber() != 0);

    //Do not parallelize this
    for (IndexType n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        if(m_mesh.isNodeFree(n))
        {
           Frho(n) = m_mesh.getNodeState(n, m_mesh.getDim() + 1);

           invMrho.diagonal()[n] = 1;
        }
    }
}

void SolverCompressible::applyBoundaryConditionsMom(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F, const Eigen::VectorXd& qVPrev)
{
    assert(m_mesh.getNodesNumber() != 0);

    const unsigned short dim = m_mesh.getDim();

    //Do not parallelize this
    for (IndexType n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
    {
        if(m_mesh.isNodeFree(n) && !m_mesh.isNodeBound(n))
        {
            for(unsigned short d = 0 ; d < dim - 1 ; ++d)
            {
                F(n + d*m_mesh.getNodesNumber()) = 0;
                invM.diagonal()[n + d*m_mesh.getNodesNumber()] = 1;
            }

             F(n + (dim - 1)*m_mesh.getNodesNumber()) = - m_gravity;
             invM.diagonal()[n + (dim - 1)*m_mesh.getNodesNumber()] = 1;
        }
        else if(m_mesh.isNodeBound(n))
        {
            std::vector<double> result;
            result = m_lua[m_mesh.getNodeType(n)](m_mesh.getNodePosition(n),
                                                  m_mesh.getNodeInitialPosition(n),
                                                  m_currentTime + m_currentDT).get<std::vector<double>>();

            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                if(m_mesh.isNodeDirichlet(n)) //Dirichlet node, directly set speed
                    F(n + d*m_mesh.getNodesNumber()) = (result[d] - qVPrev[n + d*m_mesh.getNodesNumber()])/m_currentDT;
                else                          //Not Dirichlet node, set speed through dx
                    F(n + d*m_mesh.getNodesNumber()) = ((result[d] - m_mesh.getNodePosition(n, d))/m_currentDT - qVPrev[n + d*m_mesh.getNodesNumber()])/m_currentDT;

                invM.diagonal()[n + d*m_mesh.getNodesNumber()] = 1;
            }
        }
    }
}

void SolverCompressible::displaySolverParams() const
{
    std::cout << "Initial nodes number: " << m_mesh.getNodesNumber() << std::endl;
    std::cout << "Initial elements number: " << m_mesh.getElementsNumber() << std::endl;
    std::cout << "Mesh dimension: " << m_mesh.getDim() << "D" << std::endl;
    std::cout << "alpha: " << m_mesh.getAlpha() << std::endl;
    std::cout << "hchar: " << m_mesh.getHchar() << std::endl;
    std::cout << "gamma: " << m_mesh.getGamma() << std::endl;
    std::cout << "omega: " << m_mesh.getOmega() << std::endl;
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
    std::cout << "Adapt time step: " << (m_adaptDT ? "yes" : "no") << std::endl;
    if(m_adaptDT)
    {
        std::cout << "Maximum time step: " << m_maxDT << " s" << std::endl;
        std::cout << "Security coefficient: " << m_securityCoeff << std::endl;
    }
    else
        std::cout << "Time step: " << m_maxDT << " s" << std::endl;
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

    for(unsigned short i = 0 ; i < m_pExtractor.size() ; ++i)
    {
        m_pExtractor[i]->update();
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

        solveCurrentTimeStep();

        for(unsigned short i = 0 ; i < m_pExtractor.size() ; ++i)
        {
            m_pExtractor[i]->update();
        }

        if(m_adaptDT)
        {
            double U = 0;
            double c = 0;
            for(IndexType n = 0 ; n < m_mesh.getNodesNumber() ; ++n)
            {
                if(m_mesh.isNodeBound(n) && m_mesh.isNodeFree(n))
                    continue;

                double localU = 0;
                for(unsigned short d = 0 ; d < m_mesh.getDim() ; ++d)
                    localU += m_mesh.getNodeState(n, d)*m_mesh.getNodeState(n, d);

                U = std::max(localU, U);

                c = std::max((m_K0 + m_K0prime *m_mesh.getNodeState(n, m_mesh.getDim()))/m_mesh.getNodeState(n, m_mesh.getDim() + 1), c);
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
    auto startTime = std::chrono::high_resolution_clock::now();

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

        Eigen::VectorXd qRho(m_mesh.getNodesNumber());
        qRho = invMrho*Frho;
        setNodesStatesfromQ(qRho, dim + 1, dim + 1);
        Eigen::VectorXd qP = getPFromRhoTaitMurnagham(qRho);

        setNodesStatesfromQ(qP, dim, dim);
        Frho.resize(0,0);
    }

    {
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> invM; //The mass matrix for momentum equation.
        Eigen::VectorXd F;                                  //The rhs of the momentum equation.

        buildMatricesMom(invM, F);
        applyBoundaryConditionsMom(invM, F, qVPrev);

        qAccPrev = invM*F;

        qVPrev = qV1half + 0.5*m_currentDT*qAccPrev;

        setNodesStatesfromQ(qVPrev, 0, dim - 1);
        setNodesStatesfromQ(qAccPrev, dim + 2, 2*dim + 1);
    }

    m_currentTime += m_currentDT;
    m_currentStep++;

    auto endTime = std::chrono::high_resolution_clock::now();
    auto ellapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    if(m_verboseOutput)
        std::cout << "Problem solved in " << static_cast<double>(ellapsedTime.count())/1000.0 << " s" << std::endl;

    startTime = std::chrono::high_resolution_clock::now();

    m_mesh.remesh();

    endTime = std::chrono::high_resolution_clock::now();
    ellapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    if(m_verboseOutput)
        std::cout << "Remeshing done in " << static_cast<double>(ellapsedTime.count())/1000.0 << " s" << std::endl;

    return true;
}

void SolverCompressible::buildFrho(Eigen::VectorXd& Frho)
{
    const unsigned short dim = m_mesh.getDim();

    Frho.resize(m_mesh.getNodesNumber()); Frho.setZero();

    std::vector<Eigen::VectorXd> Frhoe(m_mesh.getElementsNumber());

    #pragma omp parallel for default(shared)
    for(IndexType elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        Eigen::VectorXd Rho = getElementState(elm, dim + 1);

        Eigen::MatrixXd Mrhoe = m_mesh.getRefElementSize()*m_mesh.getElementDetJ(elm)*m_sumNTN;

        Frhoe[elm] = Mrhoe*Rho;
    }

    for(IndexType elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
            Frho(m_mesh.getElement(elm)[i]) += Frhoe[elm](i);
    }
}

void SolverCompressible::buildMatricesCont(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invMrho, Eigen::VectorXd& Frho)
{
    const unsigned short dim = m_mesh.getDim();

    invMrho.resize(m_mesh.getNodesNumber()); invMrho.setZero();

    std::vector<Eigen::MatrixXd> Mrhoe(m_mesh.getElementsNumber());
    std::vector<Eigen::VectorXd> Frhoe;

    if(!m_strongContinuity)
    {
        Frho.resize(m_mesh.getNodesNumber()); Frho.setZero();
        Frhoe.resize((m_mesh.getElementsNumber()));
    }

    #pragma omp parallel for default(shared)
    for(IndexType elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
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
            else
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

    for(IndexType elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
        {
            invMrho.diagonal()[m_mesh.getElement(elm)[i]] += Mrhoe[elm](i,i);

            if(!m_strongContinuity)
            {
                Frho(m_mesh.getElement(elm)[i]) += Frhoe[elm](i);
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(IndexType i = 0 ; i < invMrho.rows() ; ++i)
         invMrho.diagonal()[i] = 1/invMrho.diagonal()[i];
}

void SolverCompressible::buildMatricesMom(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F)
{
    const unsigned short dim = m_mesh.getDim();

    invM.resize(dim*m_mesh.getNodesNumber()); invM.setZero();
    std::vector<Eigen::MatrixXd> Me(m_mesh.getElementsNumber());

    F.resize(dim*m_mesh.getNodesNumber()); F.setZero();
    std::vector<Eigen::VectorXd> FTote(m_mesh.getElementsNumber());

    #pragma omp parallel for default(shared)
    for(IndexType elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        Eigen::VectorXd Rho = getElementState(elm, dim + 1);

        Eigen::VectorXd V((dim + 1)*dim);
        if(dim == 2)
        {
            V << getElementState(elm, 0),
                 getElementState(elm, 1);
        }
        else
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

        FTote[elm] = -Ke*V + De.transpose()*P + Fe;
    }

    for(IndexType elm = 0 ; elm < m_mesh.getElementsNumber() ; ++elm)
    {
        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                /********************************************************************
                                             Build M
                ********************************************************************/
                invM.diagonal()[m_mesh.getElement(elm)[i] + d*m_mesh.getNodesNumber()] += Me[elm](i + d*(dim + 1), i + d*(dim + 1));

                /************************************************************************
                                                Build f
                ************************************************************************/
                F(m_mesh.getElement(elm)[i] + d*m_mesh.getNodesNumber()) += FTote[elm](i + d*(dim + 1));
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(IndexType i = 0 ; i < dim*m_mesh.getNodesNumber() ; ++i)
        invM.diagonal()[i] = 1/invM.diagonal()[i];
}
