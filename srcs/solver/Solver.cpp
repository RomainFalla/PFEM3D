#include "Solver.hpp"

#include <iostream>
#include <type_traits>


Solver::Solver(const SolverCreateInfo& solverInfos) :
m_gravity(solverInfos.gravity),
m_strongPAtFS(solverInfos.strongPAtFS),
m_adaptDT(solverInfos.adaptDT),
m_currentDT(solverInfos.initialDT),
m_currentStep(0), m_currentTime(0),
m_endTime(solverInfos.endTime),
m_maxDT(solverInfos.maxDT),
m_mesh(solverInfos.meshInfos),
m_hasGMSHExtractor(false)
{
    m_solverType = SOLVER_TYPE::Undefined;

    m_lua.open_libraries(sol::lib::base, sol::lib::math);

    // set the desired number of OpenMP threads
#if defined(_OPENMP)
    const char* pNumThreads = std::getenv("OMP_NUM_THREADS");

    if(pNumThreads == nullptr)
        m_numOMPThreads = 1;
    else
        m_numOMPThreads = std::atoi(pNumThreads);

    omp_set_num_threads(m_numOMPThreads);
    Eigen::setNbThreads(m_numOMPThreads);
#else
     m_numOMPThreads = 1;
#endif

    if(m_currentDT <= 0)
        throw std::runtime_error("the initial time step should be strictly greater than 0!");

    if(m_endTime <= 0)
        throw std::runtime_error("the total time to simulate should be strictly greater than 0!");

    if(m_maxDT <= 0)
        throw std::runtime_error("the maximum time interval should be strictly greater than 0!");

    if(solverInfos.IBCfile.empty())
        throw std::runtime_error("no IBC lua file provided!");

    m_lua.script_file(solverInfos.IBCfile);
    m_lua["g"] = m_gravity;

    const uint8_t dim = m_mesh.getDim();

    if(dim == 2)
    {
        m_N1 = getNv(dim, 1);
        m_N2 = getNv(dim, 3);
        m_N3 = getNv(dim, 4);
        m_ldN2 = getN(dim - 1, 2);

        m_w1 = m_mesh.getGaussWeight(dim, 1);
        m_w2 = m_mesh.getGaussWeight(dim, 3);
        m_w3 = m_mesh.getGaussWeight(dim, 4);
        m_ldw2 = m_mesh.getGaussWeight(dim - 1, 2);

        m_m.resize(3);
        m_m << 1, 1, 0;

        m_bodyForces.resize(2);
        m_bodyForces << 0 , -m_gravity;
    }
    else
    {
        m_N1 = getNv(dim, 1);
        m_N2 = getNv(dim, 4);
        m_N3 = getNv(dim, 5);
        m_ldN2 = getN(dim - 1, 3);

        m_w1 = m_mesh.getGaussWeight(dim, 1);
        m_w2 = m_mesh.getGaussWeight(dim, 4);
        m_w3 = m_mesh.getGaussWeight(dim, 5);
        m_ldw2 = m_mesh.getGaussWeight(dim - 1, 3);

        m_m.resize(6);
        m_m << 1, 1, 1, 0, 0, 0;

        m_bodyForces.resize(3);
        m_bodyForces << 0 , 0, -m_gravity;
    }
}

void Solver::addPointExtractor(const std::string& outFileName, double timeBetweenWriting,
                               unsigned short stateToWrite, const std::vector<std::vector<double>>& points)
{
    m_pExtractor.push_back(std::make_unique<PointExtractor>(*this,
                                                            outFileName,
                                                            timeBetweenWriting,
                                                            stateToWrite,
                                                            points));
}

void Solver::addMinMaxExtractor(const std::string& outFileName, double timeBetweenWriting,
                                unsigned short coordinate, const std::string& minMax)
{
    m_pExtractor.push_back(std::make_unique<MinMaxExtractor>(*this,
                                                             outFileName,
                                                             timeBetweenWriting,
                                                             coordinate,
                                                             minMax));
}

void Solver::addMassExtractor(const std::string& outFileName, double timeBetweenWriting)
{
    m_pExtractor.push_back(std::make_unique<MassExtractor>(*this,
                                                           outFileName,
                                                           timeBetweenWriting));
}

void Solver::addGMSHExtractor(const std::string& outFileName, double timeBetweenWriting,
                              const std::vector<std::string>& whatToWrite, std::string writeAs)
{
    if(m_hasGMSHExtractor)
    {
        std::cerr << "There is already a GMSH extractor set" << std::endl;
        return;
    }

    m_pExtractor.push_back(std::make_unique<GMSHExtractor>(*this,
                                                           outFileName,
                                                           timeBetweenWriting,
                                                           whatToWrite,
                                                           writeAs));

    m_hasGMSHExtractor = true;
}

Eigen::MatrixXd Solver::getB(const Element& element) const noexcept
{
    Eigen::MatrixXd B;
    if(element.getNodeCount() == 3)
    {
        B.resize(3, 6); B.setZero();
        B(0,0) = B(2,3) = - element.getInvJ(0, 0) - element.getInvJ(1, 0);
        B(0,1) = B(2,4) = element.getInvJ(0, 0);
        B(0,2) = B(2,5) = element.getInvJ(1, 0);

        B(1,3) = B(2,0) = - element.getInvJ(0, 1) - element.getInvJ(1, 1);
        B(1,4) = B(2,1) = element.getInvJ(0, 1);
        B(1,5) = B(2,2) = element.getInvJ(1, 1);
    }
    else
    {
        B.resize(6, 12); B.setZero();
        B(0,0) = B(3,4) = B(4,8) = - element.getInvJ(0, 0) - element.getInvJ(1, 0) - element.getInvJ(2, 0);
        B(0,1) = B(3,5) = B(4,9) = element.getInvJ(0, 0);
        B(0,2) = B(3,6) = B(4,10) = element.getInvJ(1, 0);
        B(0,3) = B(3,7) = B(4,11) = element.getInvJ(2, 0);

        B(1,4) = B(3,0) = B(5,8) = - element.getInvJ(0, 1) - element.getInvJ(1, 1) - element.getInvJ(2, 1);
        B(1,5) = B(3,1) = B(5,9) = element.getInvJ(0, 1);
        B(1,6) = B(3,2) = B(5,10) = element.getInvJ(1, 1);
        B(1,7) = B(3,3) = B(5,11) = element.getInvJ(2, 1);

        B(2,8) = B(4,0) = B(5,4) = - element.getInvJ(0, 2) - element.getInvJ(1, 2) - element.getInvJ(2, 2);
        B(2,9) = B(4,1) = B(5,5) = element.getInvJ(0, 2);
        B(2,10) = B(4,2) = B(5,6) = element.getInvJ(1, 2);
        B(2,11) = B(4,3) = B(5,7) = element.getInvJ(2, 2);
    }

    return B;
}

Eigen::VectorXd Solver::getElementState(const Element& element, uint16_t state) const noexcept
{
    assert(state < m_statesNumber);

    Eigen::VectorXd stateVec(element.getNodeCount());

    if(m_mesh.getDim() == 2)
    {
        const Node& n0 = m_mesh.getNode(element.getNodeIndex(0));
        const Node& n1 = m_mesh.getNode(element.getNodeIndex(1));
        const Node& n2 = m_mesh.getNode(element.getNodeIndex(2));

        stateVec << n0.getState(state),
                    n1.getState(state),
                    n2.getState(state);
    }
    else
    {
        const Node& n0 = m_mesh.getNode(element.getNodeIndex(0));
        const Node& n1 = m_mesh.getNode(element.getNodeIndex(1));
        const Node& n2 = m_mesh.getNode(element.getNodeIndex(2));
        const Node& n3 = m_mesh.getNode(element.getNodeIndex(3));

        stateVec << n0.getState(state),
                    n1.getState(state),
                    n2.getState(state),
                    n3.getState(state);
    }

    return stateVec;
}

std::vector<Eigen::MatrixXd> Solver::getN(uint8_t dim, uint8_t nGP) const
{
    std::vector<Eigen::MatrixXd> Ns;

    std::vector<std::array<double, 3>> gP = m_mesh.getGaussPoints(dim, nGP);

    for(unsigned short p = 0 ; p < nGP ; ++p)
    {
        switch(dim)
        {
            case 1:
            {
                Eigen::MatrixXd N(1, 2); N.setZero();
                N(0,0) = (1 - gP[p][0])/2;
                N(0,1) = (1 + gP[p][0])/2;

                Ns.push_back(std::move(N));
                break;
            }

            case 2:
            {
                Eigen::MatrixXd N(1, 3); N.setZero();
                N(0,0) =  1 - gP[p][0] - gP[p][1];
                N(0,1) =  gP[p][0];
                N(0,2) =  gP[p][1];

                Ns.push_back(std::move(N));
                break;
            }

            case 3:
            {
                Eigen::MatrixXd N(1, 4); N.setZero();

                N(0,0) =  1 - gP[p][0] - gP[p][1] - gP[p][2];
                N(0,1) =  gP[p][0];
                N(0,2) =  gP[p][1];
                N(0,3) =  gP[p][2];

                Ns.push_back(std::move(N));
                break;
            }

            default:
                throw std::runtime_error("Unsupported dimension for shape functions!");
        }
    }

    return Ns;
}

std::vector<Eigen::MatrixXd> Solver::getNv(uint8_t dim, uint8_t nGP) const
{
    std::vector<Eigen::MatrixXd> Ns;

    std::vector<std::array<double, 3>> gP = m_mesh.getGaussPoints(dim, nGP);

    for(unsigned short p = 0 ; p < nGP ; ++p)
    {
        switch(dim)
        {
            case 1:
            {
                Eigen::MatrixXd N(1, 2); N.setZero();
                N(0,0) = (1 - gP[p][0])/2;
                N(0,1) = (1 + gP[p][0])/2;

                Ns.push_back(std::move(N));
                break;
            }

            case 2:
            {
                Eigen::MatrixXd N(2, 6); N.setZero();
                N(0,0) = N(1,3) =  1 - gP[p][0] - gP[p][1];
                N(0,1) = N(1,4) =  gP[p][0];
                N(0,2) = N(1,5) =  gP[p][1];

                Ns.push_back(std::move(N));
                break;
            }

            case 3:
            {
                Eigen::MatrixXd N(3, 12); N.setZero();

                N(0,0) = N(1,4) = N(2,8)  =  1 - gP[p][0] - gP[p][1] - gP[p][2];
                N(0,1) = N(1,5) = N(2,9)  =  gP[p][0];
                N(0,2) = N(1,6) = N(2,10) =  gP[p][1];
                N(0,3) = N(1,7) = N(2,11) =  gP[p][2];

                Ns.push_back(std::move(N));
                break;
            }

            default:
                throw std::runtime_error("Unsupported dimension for shape functions!");
        }
    }

    return Ns;
}

void Solver::setInitialCondition()
{
    assert(m_mesh.getNodesCount() != 0);

    //With lua, not thread safe !
    for(std::size_t n = 0 ; n < m_mesh.getNodesCount() ; ++n)
    {
        const Node& node = m_mesh.getNode(n);

        std::vector<double> result;
        result = m_lua["initStates"](node.getPosition()).get<std::vector<double>>();

        if(result.size() != m_statesNumber)
            throw std::runtime_error("Your initial condition does not set the right number of state!");

        if(node.isBound())
        {
            const auto& initFunc =  m_lua["init" + m_mesh.getNodeType(n) + "States"];
            if(initFunc.valid())
                result = m_lua["init" + m_mesh.getNodeType(n) + "States"](node.getPosition()).get<std::vector<double>>();

            bool isFixed = m_lua[m_mesh.getNodeType(n) + std::string("Fixed")].get<bool>();
            m_mesh.setNodeIsFixed(n, isFixed);
        }

        for(unsigned short i = 0 ; i < m_statesNumber ; ++i)
            m_mesh.setNodeState(n, i, result[i]);
    }
}
