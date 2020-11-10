#include "PointExtractor.hpp"

#include <Eigen/Dense>

#include "../Solver.hpp"


PointExtractor::PointExtractor(const Solver& solver, const std::string& outFileName, double timeBetweenWriting,
                               unsigned short stateToWrite, const std::vector<std::vector<double>>& points) :
Extractor(solver, outFileName, timeBetweenWriting),
m_stateToWrite(stateToWrite),
m_points(points)
{
    m_outFile.open(m_outFileName);
    if(!m_outFile.is_open())
        throw std::runtime_error("cannot open file to write point extractor: " + m_outFileName);

    if(stateToWrite > m_solver.getStatesNumber())
        throw std::runtime_error("unexpected state to write.");
}

PointExtractor::~PointExtractor()
{
    m_outFile.close();
}

void PointExtractor::update()
{
    if(m_solver.getCurrentTime() >= m_nextWriteTrigger)
    {
        const Mesh& mesh = m_solver.getMesh();

        m_outFile << std::to_string(m_solver.getCurrentTime());

        for(auto& point : m_points)
        {
            std::size_t elm;
            bool isFound = findElementIndex(mesh, elm, point);

            const Element& element = mesh.getElement(elm);

            double valueToWrite;
            if(!isFound)
               valueToWrite = 0;
            else
            {
                Eigen::MatrixXd A(mesh.getDim() + 1, mesh.getDim() + 1);
                Eigen::VectorXd b(mesh.getDim() + 1);
                if(mesh.getDim() == 2)
                {
                    const Node& n0 = mesh.getNode(element.getNodeIndex(0));
                    const Node& n1 = mesh.getNode(element.getNodeIndex(1));
                    const Node& n2 = mesh.getNode(element.getNodeIndex(2));

                    A << n0.getPosition(0), n0.getPosition(1), 1,
                         n1.getPosition(0), n1.getPosition(1), 1,
                         n2.getPosition(0), n2.getPosition(1), 1;

                    if(m_stateToWrite < m_solver.getStatesNumber())
                    {
                        b << n0.getState(m_stateToWrite),
                             n1.getState(m_stateToWrite),
                             n2.getState(m_stateToWrite);
                    }
                    else
                    {
                        std::array<double, 3> ke;
                        ke[0] = 0.5*(n0.getState(0)*n0.getState(0)
                                   + n0.getState(1)*n0.getState(1));

                        ke[1] = 0.5*(n1.getState(0)*n1.getState(0)
                                   + n1.getState(1)*n1.getState(1));

                        ke[2] = 0.5*(n2.getState(0)*n2.getState(0)
                                   + n2.getState(1)*n2.getState(1));

                        b << ke[0], ke[1], ke[2];
                    }
                }
                else
                {
                    const Node& n0 = mesh.getNode(element.getNodeIndex(0));
                    const Node& n1 = mesh.getNode(element.getNodeIndex(1));
                    const Node& n2 = mesh.getNode(element.getNodeIndex(2));
                    const Node& n3 = mesh.getNode(element.getNodeIndex(3));

                    A << n0.getPosition(0), n0.getPosition(1), n0.getPosition(2), 1,
                         n1.getPosition(0), n1.getPosition(1), n1.getPosition(2), 1,
                         n2.getPosition(0), n2.getPosition(1), n2.getPosition(2), 1,
                         n3.getPosition(0), n3.getPosition(1), n3.getPosition(2), 1;

                    if(m_stateToWrite < m_solver.getStatesNumber())
                    {
                        b << n0.getState(m_stateToWrite),
                             n1.getState(m_stateToWrite),
                             n2.getState(m_stateToWrite),
                             n3.getState(m_stateToWrite);
                    }
                    else
                    {
                        std::array<double, 4> ke;
                        ke[0] = 0.5*(n0.getState(0)*n0.getState(0)
                                   + n0.getState(1)*n0.getState(1)
                                   + n0.getState(2)*n0.getState(2));

                        ke[1] = 0.5*(n1.getState(0)*n1.getState(0)
                                   + n1.getState(1)*n1.getState(1)
                                   + n1.getState(2)*n1.getState(2));

                        ke[2] = 0.5*(n2.getState(0)*n2.getState(0)
                                   + n2.getState(1)*n2.getState(1)
                                   + n2.getState(2)*n2.getState(2));

                        ke[3] = 0.5*(n3.getState(0)*n3.getState(0)
                                   + n3.getState(1)*n3.getState(1)
                                   + n3.getState(2)*n3.getState(2));

                        b << ke[0], ke[1], ke[2], ke[3];
                    }
                }

                Eigen::VectorXd sol = A.completeOrthogonalDecomposition().solve(b);

                valueToWrite = 0;
                for(unsigned short d = 0 ; d < mesh.getDim() ; ++d)
                    valueToWrite += sol(d)*point[d];

                valueToWrite += sol(mesh.getDim());
            }
            if(mesh.getDim() == 2)
                m_outFile << "," << std::to_string(point[0]) << "," << std::to_string(point[1]) << "," << std::to_string(valueToWrite);
            else
                m_outFile << "," << std::to_string(point[0]) << "," << std::to_string(point[1]) << "," << std::to_string(point[2]) << "," << std::to_string(valueToWrite);
        }

        m_outFile << std::endl;

        m_nextWriteTrigger += m_timeBetweenWriting;
    }
}

bool PointExtractor::findElementIndex(const Mesh& mesh, std::size_t& elementIndex, const std::vector<double>& point)
{
    for(std::size_t elm = 0 ; elm < mesh.getElementsCount() ; ++elm)
    {
        const Element& element = mesh.getElement(elm);
        std::vector<double> baryCenter(mesh.getDim(), 0);

        for(uint8_t n = 0 ; n < element.getNodeCount() ; ++n)
        {
            const Node& node = mesh.getNode(element.getNodeIndex(n));
            for(uint8_t d = 0 ; d < mesh.getDim() ; ++d)
            {
                baryCenter[d] += node.getPosition(d);
            }
        }

        for(uint8_t d = 0 ; d < mesh.getDim() ; ++d)
        {
            baryCenter[d] /= static_cast<double>(element.getNodeCount());
        }

        unsigned short ok = 0;

        std::vector<double> normal(mesh.getDim());
        std::vector<double> vMiddleToBary(mesh.getDim());
        std::vector<double> vPointToVertex(mesh.getDim());

        if(mesh.getDim() == 2)
        {
            for(unsigned short i = 0 ; i < mesh.getDim() + 1 ; ++i)
            {
                unsigned short j = i + 1;
                if(i == mesh.getDim())
                    j = 0;

                const Node& ni = mesh.getNode(element.getNodeIndex(i));
                const Node& nj = mesh.getNode(element.getNodeIndex(j));

                normal[0] = nj.getPosition(1) - ni.getPosition(1);
                normal[1] = - nj.getPosition(0) + ni.getPosition(0);

                vMiddleToBary[0] = baryCenter[0] - 0.5*(ni.getPosition(0) + nj.getPosition(0));
                vMiddleToBary[1] = baryCenter[1] - 0.5*(ni.getPosition(1) + nj.getPosition(1));

                if((normal[0]*vMiddleToBary[0] + normal[1]*vMiddleToBary[1]) > 0)
                {
                    normal[0] *= -1;
                    normal[1] *= -1;
                }

                vPointToVertex[0] = point[0] - ni.getPosition(0);
                vPointToVertex[1] = point[1] - ni.getPosition(1);

                if((normal[0]*vPointToVertex[0] + normal[1]*vPointToVertex[1]) > 0)
                    break;
                else
                    ok++;
            }

            if(ok == 3)
            {
                elementIndex = elm;
                return true;
            }
        }
        else
        {
            for(unsigned short i = 0 ; i < mesh.getDim() + 1 ; ++i)
            {
                unsigned short j = i + 1;
                unsigned short k = i + 2;
                if(i == mesh.getDim() - 1)
                    k = 0;
                else if(i == mesh.getDim())
                {
                    j = 0;
                    k = 1;
                }

                const Node& ni = mesh.getNode(element.getNodeIndex(i));
                const Node& nj = mesh.getNode(element.getNodeIndex(j));
                const Node& nk = mesh.getNode(element.getNodeIndex(k));

                normal[0] = ((nj.getPosition(1) - ni.getPosition(1))
                             *(nk.getPosition(2) - ni.getPosition(2)))
                          - ((nj.getPosition(2) - ni.getPosition(2))
                             *(nk.getPosition(1) - ni.getPosition(1)));

                normal[1] = ((nj.getPosition(2) - ni.getPosition(2))
                             *(nk.getPosition(0) - ni.getPosition(0)))
                          - ((nj.getPosition(0) - ni.getPosition(0))
                             *(nk.getPosition(2) - ni.getPosition(2)));

                normal[2] = ((nj.getPosition(0) - ni.getPosition(0))
                             *(nk.getPosition(1) - ni.getPosition(1)))
                          - ((nj.getPosition(1) - ni.getPosition(1))
                             *(nk.getPosition(0) - ni.getPosition(0)));

                vMiddleToBary[0] = baryCenter[0] - (1.0f/3.0f)*(ni.getPosition(0) + nj.getPosition(0) + nk.getPosition(0));
                vMiddleToBary[1] = baryCenter[1] - (1.0f/3.0f)*(ni.getPosition(1) + nj.getPosition(1) + nk.getPosition(1));
                vMiddleToBary[2] = baryCenter[2] - (1.0f/3.0f)*(ni.getPosition(2) + nj.getPosition(2) + nk.getPosition(2));

                if((normal[0]*vMiddleToBary[0] + normal[1]*vMiddleToBary[1] + normal[2]*vMiddleToBary[2]) > 0)
                {
                    normal[0] *= -1;
                    normal[1] *= -1;
                    normal[2] *= -1;
                }

                vPointToVertex[0] = point[0] - ni.getPosition(0);
                vPointToVertex[1] = point[1] - ni.getPosition(1);
                vPointToVertex[2] = point[2] - ni.getPosition(2);

                if((normal[0]*vPointToVertex[0] + normal[1]*vPointToVertex[1] + normal[2]*vPointToVertex[2]) > 0)
                    break;
                else
                    ok++;
            }

            if(ok == 4)
            {
                elementIndex = elm;
                return true;
            }
        }
    }

    return false;
}
