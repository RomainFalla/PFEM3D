#include "MassExtractor.hpp"

#include "../Solver.hpp"


MassExtractor::MassExtractor(const Solver& solver, const std::string& outFileName, double timeBetweenWriting) :
Extractor(solver, outFileName, timeBetweenWriting)
{
    m_outFile.open(m_outFileName);
    if(!m_outFile.is_open())
        throw std::runtime_error("cannot open file to write point extractor: " + m_outFileName);
}

MassExtractor::~MassExtractor()
{
    m_outFile.close();
}

void MassExtractor::update()
{
    if(m_solver.getCurrentTime() >= m_nextWriteTrigger)
    {
        const Mesh& mesh {m_solver.getMesh()};

        double valueToWrite = 0;

        if(m_solver.getSolverType() == SOLVER_TYPE::Incompressible_PSPG)
        {
            for(std::size_t elm = 0 ; elm < mesh.getElementsCount() ; ++elm)
            {
                const Element& element = mesh.getElement(elm);
                valueToWrite += mesh.getRefElementSize(mesh.getDim())*element.getDetJ();
            }
        }
        else if(m_solver.getSolverType() == SOLVER_TYPE::WeaklyCompressible)
        {
            for(std::size_t elm = 0 ; elm < mesh.getElementsCount() ; ++elm)
            {
                const Element& element = mesh.getElement(elm);

                double volume = mesh.getRefElementSize(mesh.getDim())*element.getDetJ();

                double middleRho = 0;
                for(std::size_t n = 0 ; n < element.getNodeCount() ; ++n)
                {
                    const Node& node = mesh.getNode(element.getNodeIndex(n));
                    middleRho += node.getState(mesh.getDim() + 1);
                }

                middleRho /= static_cast<double>(element.getNodeCount());

                valueToWrite += middleRho*volume;
            }
        }
        else
            throw std::runtime_error("unsupported solver type for mass extractor!");

        m_outFile << std::to_string(m_solver.getCurrentTime()) << "," << std::to_string(valueToWrite) << std::endl;

        m_nextWriteTrigger += m_timeBetweenWriting;
    }
}
