#include "PointExtractor.hpp"

#include <Eigen/Dense>

#include "../Problem.hpp"


PointExtractor::PointExtractor(Problem* pProblem, const std::string& outFileName, double timeBetweenWriting,
                               std::string whatToWrite, const std::vector<std::vector<double>>& points) :
Extractor(pProblem, outFileName, timeBetweenWriting),
m_points(points)
{
    m_outFile.open(m_outFileName);
    if(!m_outFile.is_open())
        throw std::runtime_error("cannot open file to write point extractor: " + m_outFileName);

    std::vector<std::string> writtableData = m_pProblem->getWrittableDataName();

    if(std::find(writtableData.begin(), writtableData.end(), whatToWrite) == writtableData.end())
        throw std::runtime_error("the problem named " + m_pProblem->getID() + " cannot write data named " + whatToWrite + "!");

    if(m_pProblem->getWrittableData(whatToWrite, 0).size() != 1)
        throw std::runtime_error("PointExtractor is currently only able to extract scalars!");

    m_whatToWrite = whatToWrite;
}

PointExtractor::~PointExtractor()
{
    m_outFile.close();
}

void PointExtractor::update()
{
    if(m_pProblem->getCurrentSimTime() < m_nextWriteTrigger)
        return;

    const Mesh& mesh = m_pProblem->getMesh();

    m_outFile << std::to_string(m_pProblem->getCurrentSimTime());

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

                A << n0.getCoordinate(0), n0.getCoordinate(1), 1,
                     n1.getCoordinate(0), n1.getCoordinate(1), 1,
                     n2.getCoordinate(0), n2.getCoordinate(1), 1;

                b << m_pProblem->getWrittableData(m_whatToWrite, element.getNodeIndex(0))[0],
                     m_pProblem->getWrittableData(m_whatToWrite, element.getNodeIndex(1))[0],
                     m_pProblem->getWrittableData(m_whatToWrite, element.getNodeIndex(2))[0];
            }
            else
            {
                const Node& n0 = mesh.getNode(element.getNodeIndex(0));
                const Node& n1 = mesh.getNode(element.getNodeIndex(1));
                const Node& n2 = mesh.getNode(element.getNodeIndex(2));
                const Node& n3 = mesh.getNode(element.getNodeIndex(3));

                A << n0.getCoordinate(0), n0.getCoordinate(1), n0.getCoordinate(2), 1,
                     n1.getCoordinate(0), n1.getCoordinate(1), n1.getCoordinate(2), 1,
                     n2.getCoordinate(0), n2.getCoordinate(1), n2.getCoordinate(2), 1,
                     n3.getCoordinate(0), n3.getCoordinate(1), n3.getCoordinate(2), 1;

                b << m_pProblem->getWrittableData(m_whatToWrite, element.getNodeIndex(0))[0],
                     m_pProblem->getWrittableData(m_whatToWrite, element.getNodeIndex(1))[0],
                     m_pProblem->getWrittableData(m_whatToWrite, element.getNodeIndex(2))[0],
                     m_pProblem->getWrittableData(m_whatToWrite, element.getNodeIndex(3))[0];
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

bool PointExtractor::findElementIndex(const Mesh& mesh, std::size_t& elementIndex, const std::vector<double>& point)
{
    for(std::size_t elm = 0 ; elm < mesh.getElementsCount() ; ++elm)
    {
        const Element& element = mesh.getElement(elm);
        std::vector<double> baryCenter(mesh.getDim(), 0);

        for(unsigned int n = 0 ; n < mesh.getNodesPerElm() ; ++n)
        {
            const Node& node = mesh.getNode(element.getNodeIndex(n));
            for(unsigned int d = 0 ; d < mesh.getDim() ; ++d)
            {
                baryCenter[d] += node.getCoordinate(d);
            }
        }

        for(unsigned int d = 0 ; d < mesh.getDim() ; ++d)
        {
            baryCenter[d] /= static_cast<double>(mesh.getNodesPerElm());
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

                normal[0] = nj.getCoordinate(1) - ni.getCoordinate(1);
                normal[1] = - nj.getCoordinate(0) + ni.getCoordinate(0);

                vMiddleToBary[0] = baryCenter[0] - 0.5*(ni.getCoordinate(0) + nj.getCoordinate(0));
                vMiddleToBary[1] = baryCenter[1] - 0.5*(ni.getCoordinate(1) + nj.getCoordinate(1));

                if((normal[0]*vMiddleToBary[0] + normal[1]*vMiddleToBary[1]) > 0)
                {
                    normal[0] *= -1;
                    normal[1] *= -1;
                }

                vPointToVertex[0] = point[0] - ni.getCoordinate(0);
                vPointToVertex[1] = point[1] - ni.getCoordinate(1);

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

                normal[0] = ((nj.getCoordinate(1) - ni.getCoordinate(1))
                             *(nk.getCoordinate(2) - ni.getCoordinate(2)))
                          - ((nj.getCoordinate(2) - ni.getCoordinate(2))
                             *(nk.getCoordinate(1) - ni.getCoordinate(1)));

                normal[1] = ((nj.getCoordinate(2) - ni.getCoordinate(2))
                             *(nk.getCoordinate(0) - ni.getCoordinate(0)))
                          - ((nj.getCoordinate(0) - ni.getCoordinate(0))
                             *(nk.getCoordinate(2) - ni.getCoordinate(2)));

                normal[2] = ((nj.getCoordinate(0) - ni.getCoordinate(0))
                             *(nk.getCoordinate(1) - ni.getCoordinate(1)))
                          - ((nj.getCoordinate(1) - ni.getCoordinate(1))
                             *(nk.getCoordinate(0) - ni.getCoordinate(0)));

                vMiddleToBary[0] = baryCenter[0] - (1.0f/3.0f)*(ni.getCoordinate(0) + nj.getCoordinate(0) + nk.getCoordinate(0));
                vMiddleToBary[1] = baryCenter[1] - (1.0f/3.0f)*(ni.getCoordinate(1) + nj.getCoordinate(1) + nk.getCoordinate(1));
                vMiddleToBary[2] = baryCenter[2] - (1.0f/3.0f)*(ni.getCoordinate(2) + nj.getCoordinate(2) + nk.getCoordinate(2));

                if((normal[0]*vMiddleToBary[0] + normal[1]*vMiddleToBary[1] + normal[2]*vMiddleToBary[2]) > 0)
                {
                    normal[0] *= -1;
                    normal[1] *= -1;
                    normal[2] *= -1;
                }

                vPointToVertex[0] = point[0] - ni.getCoordinate(0);
                vPointToVertex[1] = point[1] - ni.getCoordinate(1);
                vPointToVertex[2] = point[2] - ni.getCoordinate(2);

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
