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
            for(IndexType elm = 0 ; elm < mesh.getElementsNumber() ; ++elm)
            {
                valueToWrite += mesh.getRefElementSize()*mesh.getElementDetJ(elm);
            }
        }
        else if(m_solver.getSolverType() == SOLVER_TYPE::WeaklyCompressible)
        {
            for(IndexType elm = 0 ; elm < mesh.getElementsNumber() ; ++elm)
            {
                double volume = mesh.getRefElementSize()*mesh.getElementDetJ(elm);

                double middleRho = 0;
                for(IndexType n = 0 ; n < mesh.getElement(elm).size() ; ++n)
                    middleRho += mesh.getNodeState(mesh.getElement(elm)[n], mesh.getDim() + 1);

                middleRho /= static_cast<double>(mesh.getElement(elm).size());

                valueToWrite += middleRho*volume;
            }
        }
        else
            throw std::runtime_error("unsupported solver type for mass extractor!");

        m_outFile << std::to_string(m_solver.getCurrentTime()) << "," << std::to_string(valueToWrite) << std::endl;

        m_nextWriteTrigger += m_timeBetweenWriting;
    }
}
