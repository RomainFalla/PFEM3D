#include "Problem.hpp"

#include "Solver.hpp"
#include "../../extractors/Extractors.hpp"


ProbWCompNewton::ProbWCompNewton(std::string luaFilePath) :
Problem(luaFilePath)
{
    if(m_id != "WCompNewtonNoT" && m_id != "BoussinesqWC")
        throw std::runtime_error("the ProbWCompNewton does not know id " + m_id);

    if(m_id == "WCompNewtonNoT")
        m_statesNumber = 2*m_pMesh->getDim() + 2;
    else if(m_id == "BoussinesqWC")
        m_statesNumber = 2*m_pMesh->getDim() + 3;

    m_pMesh->setStatesNumber(m_statesNumber);
    addExtractors();
    setInitialCondition();
    std::cout << "Loading solver" << std::flush;
    m_pSolver = std::make_unique<SolverWCompNewton>(this, m_pMesh.get(), m_problemParams);
    std::cout << "\rLoading solver\t\t\tok" << std::endl;
}

ProbWCompNewton::~ProbWCompNewton()
{

}

void ProbWCompNewton::displayParams() const
{
    std::cout << "Simulation time: " << m_maxTime << "\n"
              << "OpenMP threads: " << m_nThreads << "\n"
              << "Verbose output: " << std::boolalpha << m_verboseOutput << std::endl;

    m_pMesh->displayToConsole();
    m_pSolver->displayParams();
}

std::vector<std::string> ProbWCompNewton::getWrittableDataName() const
{
    if(m_pMesh->getDim() == 2)
    {
        if(m_id == "WCompNewtonNoT")
            return {"u", "v", "p", "rho", "ax", "ay", "ke", "magV", "velocity"};
        else if(m_id == "BoussinesqWC")
            return {"u", "v", "p", "rho", "ax", "ay", "T", "ke", "magV", "velocity"};
    }
    else
    {   if(m_id == "WCompNewtonNoT")
            return {"u", "v", "w", "p", "rho", "ax", "ay", "az", "ke", "magV", "velocity"};
        else if(m_id == "BoussinesqWC")
            return {"u", "v", "w", "p", "rho", "ax", "ay", "az", "T", "ke", "magV", "velocity"};
    }

    return {};
}

std::vector<double> ProbWCompNewton::getWrittableData(const std::string& name, std::size_t nodeIndex) const
{
    const Node& node = m_pMesh->getNode(nodeIndex);

    bool error = false;

    if(m_id == "WCompNewtonNoT")
    {
        if(name == "u")
            return {node.getState(0)};
        else if(name == "v")
            return {node.getState(1)};
        else
        {
            if(m_pMesh->getDim() == 2)
            {
                if(name == "p")
                    return {node.getState(2)};
                else if(name == "rho")
                    return {node.getState(3)};
                else if(name == "ax")
                    return {node.getState(4)};
                else if(name == "ay")
                    return {node.getState(5)};
                else if(name == "ke")
                    return {0.5*(node.getState(0)*node.getState(0) + node.getState(1)*node.getState(1))};
                else if(name == "magV")
                    return {std::sqrt(node.getState(0)*node.getState(0) + node.getState(1)*node.getState(1))};
                else if(name == "velocity")
                    return {node.getState(0), node.getState(1), 0};
                else
                    error = true;

            }
            else
            {
                if(name == "w")
                    return {node.getState(2)};
                else if(name == "p")
                    return {node.getState(3)};
                else if(name == "rho")
                    return {node.getState(4)};
                else if(name == "ax")
                    return {node.getState(5)};
                else if(name == "ay")
                    return {node.getState(6)};
                else if(name == "ay")
                    return {node.getState(7)};
                else if(name == "ke")
                    return {0.5*(node.getState(0)*node.getState(0) + node.getState(1)*node.getState(1) + node.getState(2)*node.getState(2))};
                else if(name == "magV")
                    return {std::sqrt(node.getState(0)*node.getState(0) + node.getState(1)*node.getState(1) + node.getState(2)*node.getState(2))};
                else if(name == "velocity")
                    return {node.getState(0), node.getState(1), node.getState(2)};
                else
                    error = true;
            }
        }
    }
    else if(m_id == "BoussinesqWC")
    {
        if(name == "u")
            return {node.getState(0)};
        else if(name == "v")
            return {node.getState(1)};
        else
        {
            if(m_pMesh->getDim() == 2)
            {
                if(name == "p")
                    return {node.getState(2)};
                else if(name == "rho")
                    return {node.getState(3)};
                else if(name == "ax")
                    return {node.getState(4)};
                else if(name == "ay")
                    return {node.getState(5)};
                else if(name == "T")
                    return {node.getState(6)};
                else if(name == "ke")
                    return {0.5*(node.getState(0)*node.getState(0) + node.getState(1)*node.getState(1))};
                else if(name == "magV")
                    return {std::sqrt(node.getState(0)*node.getState(0) + node.getState(1)*node.getState(1))};
                else if(name == "velocity")
                    return {node.getState(0), node.getState(1), 0};
                else
                    error = true;

            }
            else
            {
                if(name == "w")
                    return {node.getState(2)};
                else if(name == "p")
                    return {node.getState(3)};
                else if(name == "rho")
                    return {node.getState(4)};
                else if(name == "ax")
                    return {node.getState(5)};
                else if(name == "ay")
                    return {node.getState(6)};
                else if(name == "ay")
                    return {node.getState(7)};
                else if(name == "T")
                    return {node.getState(8)};
                else if(name == "ke")
                    return {0.5*(node.getState(0)*node.getState(0) + node.getState(1)*node.getState(1) + node.getState(2)*node.getState(2))};
                else if(name == "magV")
                    return {std::sqrt(node.getState(0)*node.getState(0) + node.getState(1)*node.getState(1) + node.getState(2)*node.getState(2))};
                else if(name == "velocity")
                    return {node.getState(0), node.getState(1), node.getState(2)};
                else
                    error = true;
            }
        }
    }

    if(error)
        throw std::runtime_error("The data " + name + " cannot be extract from problem " + getID());

    return {}; //silence waring
}

std::vector<std::string> ProbWCompNewton::getGlobalWrittableDataName() const
{
    return {"mass"};
}

double ProbWCompNewton::getGlobalWrittableData(const std::string& name) const
{
    if(name == "mass")
    {
        double mass = 0;

        for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
        {
            const Element& element = m_pMesh->getElement(elm);

            double volume = element.getSize();

            double middleRho = 0;
            for(std::size_t n = 0 ; n < m_pMesh->getNodesPerElm() ; ++n)
            {
                const Node& node = m_pMesh->getNode(element.getNodeIndex(n));
                middleRho += node.getState(m_pMesh->getDim() + 1);
            }

            middleRho /= static_cast<double>(m_pMesh->getNodesPerElm());

            mass += middleRho*volume;
        }

        return mass;
    }
    else
        throw std::runtime_error("The global data " + name + " cannot be extract from problem " + getID());
}
