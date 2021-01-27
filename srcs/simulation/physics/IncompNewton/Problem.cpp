#include "Problem.hpp"

#include "Solver.hpp"
#include "../../extractors/Extractors.hpp"


ProbIncompNewton::ProbIncompNewton(std::string luaFilePath) :
Problem(luaFilePath)
{
    if(m_id != "IncompNewtonNoT" && m_id != "Boussinesq" && m_id != "Conduction")
        throw std::runtime_error("the ProbIncompNewton does not know id " + m_id);

    if(m_id == "IncompNewtonNoT")
        m_statesNumber = m_pMesh->getDim() + 1;
    else if(m_id == "Boussinesq")
        m_statesNumber = m_pMesh->getDim() + 2;
    else if(m_id == "Conduction")
        m_statesNumber = 1;

    m_pMesh->setStatesNumber(m_statesNumber);
    addExtractors();
    setInitialCondition();
    std::cout << "Loading solver" << std::flush;
    m_pSolver = std::make_unique<SolverIncompNewton>(this, m_pMesh.get(), m_problemParams);
    std::cout << "\rLoading solver\t\t\tok" << std::endl;
}

ProbIncompNewton::~ProbIncompNewton()
{

}

void ProbIncompNewton::displayParams() const
{
    std::cout << "Simulation time: " << m_maxTime << "\n"
              << "OpenMP threads: " << m_nThreads << "\n"
              << "Verbose output: " << std::boolalpha << m_verboseOutput << std::endl;

    m_pMesh->displayToConsole();
    m_pSolver->displayParams();
}

std::vector<std::string> ProbIncompNewton::getWrittableDataName() const
{
    if(m_pMesh->getDim() == 2)
    {
        if(m_id == "IncompNewtonNoT")
            return {"u", "v", "p", "ke", "magV", "velocity"};
        else if(m_id == "Boussinesq")
            return {"u", "v", "p", "T", "ke", "magV", "velocity"};
        else if(m_id == "Conduction")
            return {"T"};
    }
    else
    {   if(m_id == "IncompNewtonNoT")
            return {"u", "v", "w", "p", "ke", "magV", "velocity"};
        else if(m_id == "Boussinesq")
            return {"u", "v", "w", "p", "T", "ke", "magV", "velocity"};
        else if(m_id == "Conduction")
            return {"T"};
    }

    return {};
}

std::vector<double> ProbIncompNewton::getWrittableData(const std::string& name, std::size_t nodeIndex) const
{
    const Node& node = m_pMesh->getNode(nodeIndex);

    bool error = false;

    if(m_id == "IncompNewtonNoT")
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
    else if(m_id == "Boussinesq")
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
                else if(name == "T")
                    return {node.getState(3)};
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
                else if(name == "T")
                    return {node.getState(4)};
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
    else if(m_id == "Conduction")
    {
        if(name == "T")
            return {node.getState(0)};
        else
            error = true;
    }

    if(error)
        throw std::runtime_error("The data " + name + " cannot be extract from problem " + getID());

    return {}; //silence waring
}

std::vector<std::string> ProbIncompNewton::getGlobalWrittableDataName() const
{
    return {"mass"};
}

double ProbIncompNewton::getGlobalWrittableData(const std::string& name) const
{
    if(name == "mass")
    {
        double mass = 0;

        for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
        {
            const Element& element = m_pMesh->getElement(elm);
            mass += element.getSize();
        }

        return mass;
    }
    else
        throw std::runtime_error("The global data " + name + " cannot be extract from problem " + getID());
}
