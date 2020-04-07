#include "PointExtractor.hpp"

#include <iostream>
#include <Eigen/Dense>


PointExtractor::PointExtractor(const std::string& outFileName, double timeBetweenWriting,
                               unsigned short stateToWrite, const std::vector<std::vector<double>>& points,
                               unsigned short statesNumber) :
Extractor(outFileName, timeBetweenWriting), m_statesNumber(statesNumber), m_stateToWrite(stateToWrite), m_points(points)
{
    m_outFile.open(m_outFileName);
    if(!m_outFile.is_open())
    {
        std::string errorText = std::string("cannot open file to write point extractor: ") + m_outFileName;
        throw std::runtime_error(errorText);
    }
    if(stateToWrite > statesNumber)
        throw std::runtime_error("unexpected state to write.");
}

PointExtractor::~PointExtractor()
{
    m_outFile.close();
}

void PointExtractor::update(const Mesh& mesh, double currentTime, unsigned int currentStep)
{
    if(currentTime >= m_nextWriteTrigger)
    {
        m_outFile << std::to_string(currentTime);
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
                    A << mesh.getNodePosition(mesh.getElement(elm)[0], 0), mesh.getNodePosition(mesh.getElement(elm)[0], 1), 1,
                         mesh.getNodePosition(mesh.getElement(elm)[1], 0), mesh.getNodePosition(mesh.getElement(elm)[1], 1), 1,
                         mesh.getNodePosition(mesh.getElement(elm)[2], 0), mesh.getNodePosition(mesh.getElement(elm)[2], 1), 1;
                }
                else
                {
                    A << mesh.getNodePosition(mesh.getElement(elm)[0], 0), mesh.getNodePosition(mesh.getElement(elm)[0], 1), mesh.getNodePosition(mesh.getElement(elm)[0], 2), 1,
                         mesh.getNodePosition(mesh.getElement(elm)[1], 0), mesh.getNodePosition(mesh.getElement(elm)[1], 1), mesh.getNodePosition(mesh.getElement(elm)[1], 2), 1,
                         mesh.getNodePosition(mesh.getElement(elm)[2], 0), mesh.getNodePosition(mesh.getElement(elm)[2], 1), mesh.getNodePosition(mesh.getElement(elm)[2], 2), 1,
                         mesh.getNodePosition(mesh.getElement(elm)[3], 0), mesh.getNodePosition(mesh.getElement(elm)[3], 1), mesh.getNodePosition(mesh.getElement(elm)[3], 2), 1;
                }

                Eigen::VectorXd b(mesh.getDim() + 1);
                if(m_stateToWrite < m_statesNumber)
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

            m_outFile << "," << std::to_string(point[0]) << "," << std::to_string(point[1]) << "," << std::to_string(valueToWrite);
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

        if(mesh.getDim() == 2)
        {
            std::vector<double> normal(mesh.getDim());
            std::vector<double> vMiddleToBary(mesh.getDim());
            std::vector<double> vPointToVertex(mesh.getDim());

            for(unsigned short i = 0 ; i < mesh.getDim() + 1 ; ++i)
            {
                //Segment 1-0
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
            throw std::runtime_error("3D meshes currently unsupported for Point Extractor");

    }

    return false;
}
