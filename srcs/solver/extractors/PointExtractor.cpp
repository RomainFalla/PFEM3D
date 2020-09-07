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
            IndexType elm;
            bool isFound = findElementIndex(mesh, elm, point);

            double valueToWrite;
            if(!isFound)
               valueToWrite = 0;
            else
            {
                Eigen::MatrixXd A(mesh.getDim() + 1, mesh.getDim() + 1);
                if(mesh.getDim() == 2)
                {
                    std::vector<double> posV0 = mesh.getNodePosition(mesh.getElement(elm)[0]);
                    std::vector<double> posV1 = mesh.getNodePosition(mesh.getElement(elm)[1]);
                    std::vector<double> posV2 = mesh.getNodePosition(mesh.getElement(elm)[2]);

                    A << posV0[0], posV0[1], 1,
                         posV1[0], posV1[1], 1,
                         posV2[0], posV2[1], 1;
                }
                else
                {
                    std::vector<double> posV0 = mesh.getNodePosition(mesh.getElement(elm)[0]);
                    std::vector<double> posV1 = mesh.getNodePosition(mesh.getElement(elm)[1]);
                    std::vector<double> posV2 = mesh.getNodePosition(mesh.getElement(elm)[2]);
                    std::vector<double> posV3 = mesh.getNodePosition(mesh.getElement(elm)[3]);

                    A << posV0[0], posV0[1], posV0[2], 1,
                         posV1[0], posV1[1], posV1[2], 1,
                         posV2[0], posV2[1], posV2[2], 1,
                         posV3[0], posV3[1], posV3[2], 1;
                }

                Eigen::VectorXd b(mesh.getDim() + 1);
                if(m_stateToWrite < m_solver.getStatesNumber())
                {
                    if(mesh.getDim() == 2)
                    {
                        b << mesh.getNodeState(mesh.getElement(elm)[0], m_stateToWrite),
                             mesh.getNodeState(mesh.getElement(elm)[1], m_stateToWrite),
                             mesh.getNodeState(mesh.getElement(elm)[2], m_stateToWrite);
                    }
                    else
                    {
                        b << mesh.getNodeState(mesh.getElement(elm)[0], m_stateToWrite),
                             mesh.getNodeState(mesh.getElement(elm)[1], m_stateToWrite),
                             mesh.getNodeState(mesh.getElement(elm)[2], m_stateToWrite),
                             mesh.getNodeState(mesh.getElement(elm)[3], m_stateToWrite);
                    }
                }
                else
                {
                    std::vector<double> ke(mesh.getDim() + 1);
                    for(unsigned short i = 0 ; i < ke.size() ; ++i)
                    {
                        ke[i] = 0;
                        for(unsigned short d = 0 ; d < mesh.getDim() ; ++d)
                        {
                            ke[i]+= mesh.getNodeState(mesh.getElement(elm)[i], d)*
                                    mesh.getNodeState(mesh.getElement(elm)[i], d);
                        }
                        ke[i] = 0.5*std::sqrt(ke[i]);
                    }

                    if(mesh.getDim() == 2)
                    {
                        b << ke[0],
                             ke[1],
                             ke[2];
                    }
                    else
                    {
                         b << ke[0],
                              ke[1],
                              ke[2],
                              ke[3];
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

bool PointExtractor::findElementIndex(const Mesh& mesh, IndexType& elementIndex, const std::vector<double>& point)
{
    for(IndexType elm = 0 ; elm < mesh.getElementsNumber() ; ++elm)
    {
        std::vector<double> baryCenter(mesh.getDim());
        for(unsigned short d = 0 ; d < mesh.getDim() ; ++d)
        {
            baryCenter[d] = 0;
            for(unsigned short n = 0 ; n < mesh.getDim() + 1 ; ++n)
            {
                baryCenter[d] += mesh.getNodePosition(mesh.getElement(elm)[n], d);
            }
            baryCenter[d] /= static_cast<double>(mesh.getDim() + 1);
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

                normal[0] = mesh.getNodePosition(mesh.getElement(elm)[j], 1) - mesh.getNodePosition(mesh.getElement(elm)[i], 1);
                normal[1] = - mesh.getNodePosition(mesh.getElement(elm)[j], 0) + mesh.getNodePosition(mesh.getElement(elm)[i], 0);

                vMiddleToBary[0] = baryCenter[0] - 0.5*(mesh.getNodePosition(mesh.getElement(elm)[i], 0) + mesh.getNodePosition(mesh.getElement(elm)[j], 0));
                vMiddleToBary[1] = baryCenter[1] - 0.5*(mesh.getNodePosition(mesh.getElement(elm)[i], 1) + mesh.getNodePosition(mesh.getElement(elm)[j], 1));

                if((normal[0]*vMiddleToBary[0] + normal[1]*vMiddleToBary[1]) > 0)
                {
                    normal[0] *= -1;
                    normal[1] *= -1;
                }

                vPointToVertex[0] = point[0] - mesh.getNodePosition(mesh.getElement(elm)[i], 0);
                vPointToVertex[1] = point[1] - mesh.getNodePosition(mesh.getElement(elm)[i], 1);

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

                normal[0] = ((mesh.getNodePosition(mesh.getElement(elm)[j], 1) - mesh.getNodePosition(mesh.getElement(elm)[i], 1))
                             *(mesh.getNodePosition(mesh.getElement(elm)[k], 2) - mesh.getNodePosition(mesh.getElement(elm)[i], 2)))
                          - ((mesh.getNodePosition(mesh.getElement(elm)[j], 2) - mesh.getNodePosition(mesh.getElement(elm)[i], 2))
                             *(mesh.getNodePosition(mesh.getElement(elm)[k], 1) - mesh.getNodePosition(mesh.getElement(elm)[i], 1))) ;

                normal[1] = ((mesh.getNodePosition(mesh.getElement(elm)[j], 2) - mesh.getNodePosition(mesh.getElement(elm)[i], 2))
                             *(mesh.getNodePosition(mesh.getElement(elm)[k], 0) - mesh.getNodePosition(mesh.getElement(elm)[i], 0)))
                          - ((mesh.getNodePosition(mesh.getElement(elm)[j], 0) - mesh.getNodePosition(mesh.getElement(elm)[i], 0))
                             *(mesh.getNodePosition(mesh.getElement(elm)[k], 2) - mesh.getNodePosition(mesh.getElement(elm)[i], 2))) ;

                normal[2] = ((mesh.getNodePosition(mesh.getElement(elm)[j], 0) - mesh.getNodePosition(mesh.getElement(elm)[i], 0))
                             *(mesh.getNodePosition(mesh.getElement(elm)[k], 1) - mesh.getNodePosition(mesh.getElement(elm)[i], 1)))
                          - ((mesh.getNodePosition(mesh.getElement(elm)[j], 1) - mesh.getNodePosition(mesh.getElement(elm)[i], 1))
                             *(mesh.getNodePosition(mesh.getElement(elm)[k], 0) - mesh.getNodePosition(mesh.getElement(elm)[i], 0))) ;

                vMiddleToBary[0] = baryCenter[0] - (1.0f/3.0f)*(mesh.getNodePosition(mesh.getElement(elm)[i], 0) + mesh.getNodePosition(mesh.getElement(elm)[j], 0) + mesh.getNodePosition(mesh.getElement(elm)[k], 0));
                vMiddleToBary[1] = baryCenter[1] - (1.0f/3.0f)*(mesh.getNodePosition(mesh.getElement(elm)[i], 1) + mesh.getNodePosition(mesh.getElement(elm)[j], 1) + mesh.getNodePosition(mesh.getElement(elm)[k], 1));
                vMiddleToBary[2] = baryCenter[2] - (1.0f/3.0f)*(mesh.getNodePosition(mesh.getElement(elm)[i], 2) + mesh.getNodePosition(mesh.getElement(elm)[j], 2) + mesh.getNodePosition(mesh.getElement(elm)[k], 2));

                if((normal[0]*vMiddleToBary[0] + normal[1]*vMiddleToBary[1] + normal[2]*vMiddleToBary[2]) > 0)
                {
                    normal[0] *= -1;
                    normal[1] *= -1;
                    normal[2] *= -1;
                }

                vPointToVertex[0] = point[0] - mesh.getNodePosition(mesh.getElement(elm)[i], 0);
                vPointToVertex[1] = point[1] - mesh.getNodePosition(mesh.getElement(elm)[i], 1);
                vPointToVertex[2] = point[2] - mesh.getNodePosition(mesh.getElement(elm)[i], 2);

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
