#include "Mesh.hpp"

#include <fstream>
#include <iostream>

#include <gmsh.h>


Mesh::Mesh(const nlohmann::json& j)
{
    m_verboseOutput = j["verboseOutput"].get<bool>();

    m_hchar       = j["Remeshing"]["hchar"].get<double>();
    m_alpha       = j["Remeshing"]["alpha"].get<double>();
    m_omega       = j["Remeshing"]["omega"].get<double>();
    m_gamma       = j["Remeshing"]["gamma"].get<double>();
    m_boundingBox = j["Remeshing"]["boundingBox"].get<std::vector<double>>();
}

Mesh::~Mesh()
{

}

bool Mesh::addNodes()
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh!");
    assert(m_elementsDetJ.size() == m_elementsList.size() && "The determinant are not computed!");

    bool addedNodes = false;

    for(std::size_t i = 0 ; i < m_elementsList.size() ; ++i)
    {
        //If an element is too big, we add a node at his center
        if(m_elementsDetJ[i]*getRefElementSize() > m_omega*std::pow(m_hchar, m_dim))
        {
            Node newNode(m_dim);

            for(unsigned short k = 0 ; k < m_dim ; ++k)
            {
                newNode.position[k] = 0;
                for(unsigned short d = 0 ; d <= m_dim ; ++d)
                {
                    newNode.position[k] += m_nodesList[m_elementsList[i][d]].position[k];
                }
                newNode.position[k] /= (m_dim + 1);
            }

            newNode.states.resize(m_nodesList[0].states.size());
            for(unsigned short k = 0 ; k < m_nodesList[0].states.size() ; ++k)
            {
                newNode.states[k] = 0;
                for(unsigned short d = 0 ; d <= m_dim ; ++d)
                {
                    newNode.states[k] += m_nodesList[m_elementsList[i][d]].states[k];
                }
                newNode.states[k] /= (m_dim + 1);
            }

            m_nodesList.push_back(newNode);

            if(m_verboseOutput)
            {
                std::cout << "Adding node (" << "(";
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                {
                    std::cout << newNode.position[d];
                    if(d == m_dim - 1)
                        std::cout << ")";
                    else
                        std::cout << ", ";
                }
                std::cout << std::endl;
            }

            addedNodes = true;
        }
    }

    return addedNodes;
}

bool Mesh::checkBoundingBox()
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh !");

    bool outofBBNodes = false;

    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        //If the node is out of the bounding box, we delete it.
        // Bounding box fromat: [xmin, ymin, zmin, xmax, ymax, zmax]
        for(unsigned short d = 0 ; d < m_dim ; ++d)
        {
            if(m_nodesList[n].position[d] < m_boundingBox[d] ||
               m_nodesList[n].position[d] > m_boundingBox[d + m_dim])
            {
                m_nodesList[n].toBeDeleted = true;
                outofBBNodes = true;
                break;
            }
        }
    }

    m_nodesList.erase(
    std::remove_if(m_nodesList.begin(), m_nodesList.end(), [this](const Node& node)
    {
       if(node.toBeDeleted)
       {
           if(this->m_verboseOutput)
           {
               std::cout << "Removing out of bounding box node (";
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                {
                    std::cout << node.position[d];
                    if(d == m_dim - 1)
                        std::cout << ")";
                    else
                        std::cout << ", ";
                }
                std::cout << std::endl;
           }

           return true;
       }
       else
           return false;
    }), m_nodesList.end());

    return outofBBNodes;
}

void Mesh::computeElementsDetJ()
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh !");
    assert(!m_elementsJ.empty() && m_elementsJ.size() == m_elementsList.size());

    m_elementsDetJ.clear();
    m_elementsDetJ.resize(m_elementsList.size());

    #pragma omp parallel for default(shared)
    for(std::size_t i = 0 ; i < m_elementsList.size() ; ++i)
    {
        if(m_dim == 2)
        {
            m_elementsDetJ[i] = m_elementsJ[i][0][0]*m_elementsJ[i][1][1]
                              - m_elementsJ[i][1][0]*m_elementsJ[i][0][1];
        }
        else if(m_dim == 3)
        {
            m_elementsDetJ[i] = m_elementsJ[i][0][0]*m_elementsJ[i][1][1]*m_elementsJ[i][2][2]
                              + m_elementsJ[i][0][1]*m_elementsJ[i][1][2]*m_elementsJ[i][2][0]
                              + m_elementsJ[i][0][2]*m_elementsJ[i][1][0]*m_elementsJ[i][2][1]
                              - m_elementsJ[i][2][0]*m_elementsJ[i][1][1]*m_elementsJ[i][0][2]
                              - m_elementsJ[i][2][1]*m_elementsJ[i][1][2]*m_elementsJ[i][0][0]
                              - m_elementsJ[i][2][2]*m_elementsJ[i][1][0]*m_elementsJ[i][0][1];
        }
    }
}

void Mesh::computeElementsInvJ()
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh!");
    assert(!m_elementsJ.empty() && m_elementsJ.size() == m_elementsList.size());
    assert(m_elementsDetJ.size() == m_elementsList.size() && "The determinant are not computed!");

    m_elementsInvJ.clear();
    m_elementsInvJ.resize(m_elementsList.size());

    #pragma omp parallel for default(shared)
    for(std::size_t i = 0 ; i < m_elementsList.size() ; ++i)
    {
        if(m_dim == 2)
        {
            m_elementsInvJ[i].resize(2);
            m_elementsInvJ[i][0].resize(2);
            m_elementsInvJ[i][1].resize(2);

            m_elementsInvJ[i][0][0] = m_elementsJ[i][1][1]/m_elementsDetJ[i];

            m_elementsInvJ[i][0][1] = - m_elementsJ[i][0][1]/m_elementsDetJ[i];

            m_elementsInvJ[i][1][0] = - m_elementsJ[i][1][0]/m_elementsDetJ[i];

            m_elementsInvJ[i][1][1] = m_elementsJ[i][0][0]/m_elementsDetJ[i];
        }
        else if(m_dim == 3)
        {
            m_elementsInvJ[i].resize(3);
            m_elementsInvJ[i][0].resize(3);
            m_elementsInvJ[i][1].resize(3);
            m_elementsInvJ[i][2].resize(3);

            m_elementsInvJ[i][0][0] = (m_elementsJ[i][1][1]*m_elementsJ[i][2][2]
                                    - m_elementsJ[i][1][2]*m_elementsJ[i][2][1])/m_elementsDetJ[i];

            m_elementsInvJ[i][0][1] = (m_elementsJ[i][2][1]*m_elementsJ[i][0][2]
                                    - m_elementsJ[i][2][2]*m_elementsJ[i][0][1])/m_elementsDetJ[i];

            m_elementsInvJ[i][0][2] = (m_elementsJ[i][0][1]*m_elementsJ[i][1][2]
                                    - m_elementsJ[i][0][2]*m_elementsJ[i][1][1])/m_elementsDetJ[i];

            m_elementsInvJ[i][1][0] = (m_elementsJ[i][2][0]*m_elementsJ[i][1][2]
                                    - m_elementsJ[i][1][0]*m_elementsJ[i][2][2])/m_elementsDetJ[i];

            m_elementsInvJ[i][1][1] = (m_elementsJ[i][0][0]*m_elementsJ[i][2][2]
                                    - m_elementsJ[i][2][0]*m_elementsJ[i][0][2])/m_elementsDetJ[i];

            m_elementsInvJ[i][1][2] = (m_elementsJ[i][1][0]*m_elementsJ[i][0][2]
                                    - m_elementsJ[i][0][0]*m_elementsJ[i][1][2])/m_elementsDetJ[i];

            m_elementsInvJ[i][2][0] = (m_elementsJ[i][1][0]*m_elementsJ[i][2][1]
                                    - m_elementsJ[i][2][0]*m_elementsJ[i][1][1])/m_elementsDetJ[i];

            m_elementsInvJ[i][2][1] = (m_elementsJ[i][2][0]*m_elementsJ[i][0][1]
                                    - m_elementsJ[i][0][0]*m_elementsJ[i][2][1])/m_elementsDetJ[i];

            m_elementsInvJ[i][2][2] = (m_elementsJ[i][0][0]*m_elementsJ[i][1][1]
                                    - m_elementsJ[i][1][0]*m_elementsJ[i][0][1])/m_elementsDetJ[i];
        }

    }
}

void Mesh::computeElementsJ()
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh!");

    m_elementsJ.clear();
    m_elementsJ.resize(m_elementsList.size());

    #pragma omp parallel for default(shared)
    for(std::size_t i = 0 ; i < m_elementsList.size() ; ++i)
    {
        if(m_dim == 2)
        {
            m_elementsJ[i].resize(2);
            m_elementsJ[i][0].resize(2);
            m_elementsJ[i][1].resize(2);

            double x0 = m_nodesList[m_elementsList[i][0]].position[0];
            double x1 = m_nodesList[m_elementsList[i][1]].position[0];
            double x2 = m_nodesList[m_elementsList[i][2]].position[0];
            double y0 = m_nodesList[m_elementsList[i][0]].position[1];
            double y1 = m_nodesList[m_elementsList[i][1]].position[1];
            double y2 = m_nodesList[m_elementsList[i][2]].position[1];

            m_elementsJ[i][0][0] = x1 - x0;
            m_elementsJ[i][0][1] = x2 - x0;
            m_elementsJ[i][1][0] = y1 - y0;
            m_elementsJ[i][1][1] = y2 - y0;
        }
        else if(m_dim == 3)
        {
            m_elementsJ[i].resize(3);
            m_elementsJ[i][0].resize(3);
            m_elementsJ[i][1].resize(3);
            m_elementsJ[i][2].resize(3);

            double x0 = m_nodesList[m_elementsList[i][0]].position[0];
            double x1 = m_nodesList[m_elementsList[i][1]].position[0];
            double x2 = m_nodesList[m_elementsList[i][2]].position[0];
            double x3 = m_nodesList[m_elementsList[i][3]].position[0];
            double y0 = m_nodesList[m_elementsList[i][0]].position[1];
            double y1 = m_nodesList[m_elementsList[i][1]].position[1];
            double y2 = m_nodesList[m_elementsList[i][2]].position[1];
            double y3 = m_nodesList[m_elementsList[i][3]].position[1];
            double z0 = m_nodesList[m_elementsList[i][0]].position[2];
            double z1 = m_nodesList[m_elementsList[i][1]].position[2];
            double z2 = m_nodesList[m_elementsList[i][2]].position[2];
            double z3 = m_nodesList[m_elementsList[i][3]].position[2];

            m_elementsJ[i][0][0] = x1 - x0;
            m_elementsJ[i][0][1] = x2 - x0;
            m_elementsJ[i][0][2] = x3 - x0;
            m_elementsJ[i][1][0] = y1 - y0;
            m_elementsJ[i][1][1] = y2 - y0;
            m_elementsJ[i][1][2] = y3 - y0;
            m_elementsJ[i][2][0] = z1 - z0;
            m_elementsJ[i][2][1] = z2 - z0;
            m_elementsJ[i][2][2] = z3 - z0;
        }
    }
}

//void Mesh::computeFreeSurfaceEdgeDetJ()
//{
//    assert(!m_freeSurfaceEdgesList.empty() && !m_nodesList.empty() && "There is no mesh !");
//
//    m_freeSurfaceEdgeDetJ.clear();
//    m_freeSurfaceEdgeDetJ.resize(m_freeSurfaceEdgesList.size());
//
//    #pragma omp parallel for default(shared)
//    for(std::size_t i = 0 ; i < m_freeSurfaceEdgesList.size() ; ++i)
//    {
//        m_freeSurfaceEdgeDetJ[i] = 0; //TO DO
//    }
//}

unsigned short Mesh::computeMeshDim() const
{
    int elementDim = -1;

    // loop over the dimension i to get the maximum element dimension in the mesh
    for(unsigned short i = 0 ; i <= 3 ; ++i)
    {
        std::vector<int> eleTypes;
        gmsh::model::mesh::getElementTypes(eleTypes, i);

        switch(eleTypes.size())
        {
            case 0:
                break;
            case 1:
                elementDim = i;
                break;
            default:
                elementDim = i;
                std::cerr   << "Hybrid meshes not handled in this example!"
                            << std::endl;
        }
    }

    return static_cast<unsigned short>(elementDim);
}

void Mesh::loadFromFile(const std::string& fileName)
{
    m_nodesList.clear();

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);

    std::ifstream file(fileName);
    if(file.is_open())
        file.close();
    else
        throw std::runtime_error("The input .msh file does not exist!");

    gmsh::open(fileName);

    // Check that the mesh is not 3D
    m_dim = computeMeshDim();

    if(m_boundingBox.size() != 2*m_dim)
        throw std::runtime_error("Bad bounding box format! Format: [xmin, ymin, xmax, ymax]");

    // We retrieve the tags of the physical groups of dimension m_dim and
    // m_dim-1
    std::vector<std::pair<int, int>> physGroupHandlesHD;
    gmsh::model::getPhysicalGroups(physGroupHandlesHD, m_dim);

    std::vector<std::pair<int, int>> physGroupHandlesLD;
    gmsh::model::getPhysicalGroups(physGroupHandlesLD, m_dim - 1);

    //The file should contain physical group for the boundary, the free surface and the fluid

    for(auto physGroupLD : physGroupHandlesLD)
    {
        std::string name;
        gmsh::model::getPhysicalName(m_dim - 1, physGroupLD.second, name);
        if(name == "Boundary" || name == "FluidInput")
        {
            std::vector<double> coord;
            std::vector<std::size_t> dummyNodesTagsBoundary;
            gmsh::model::mesh::getNodesForPhysicalGroup(m_dim - 1, physGroupLD.second,
                                                        dummyNodesTagsBoundary, coord);

            for(std::size_t i = 0 ; i < dummyNodesTagsBoundary.size() ; ++i)
            {
                Node node(m_dim);
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                    node.position[d] = coord[3*i + d];

                //If the nodes is already on the boundary, we do not add it twice
                if(std::find(m_nodesList.begin(), m_nodesList.end(), node) != m_nodesList.end())
                    continue;

                node.isBound = true;
                if(name == "FluidInput")
                    node.isFluidInput = true;

                m_nodesList.push_back(node);

                if(m_verboseOutput)
                {
                    std::cout << "Loading boundary node: " << "(";
                    for(unsigned short d = 0 ; d < m_dim ; ++d)
                    {
                        std::cout << node.position[d];
                        if(d == m_dim - 1)
                            std::cout << ")";
                        else
                            std::cout << ", ";
                    }
                    std::cout << std::endl;
                }
            }
        }
    }

    for(auto physGroupHD : physGroupHandlesHD)
    {
        std::string name;
        gmsh::model::getPhysicalName(m_dim, physGroupHD.second, name);
        if(name == "Fluid")
        {
            std::vector<std::size_t> dummyNodesTags;
            std::vector<double> coord;
            gmsh::model::mesh::getNodesForPhysicalGroup(m_dim, physGroupHD.second,
                                                        dummyNodesTags, coord);

            for(std::size_t i = 0 ; i < dummyNodesTags.size() ; ++i)
            {
                Node node(m_dim);
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                    node.position[d] = coord[3*i + d];

                //If the nodes is already on the boundary, we do not add it twice
                if(std::find(m_nodesList.begin(), m_nodesList.end(), node) != m_nodesList.end())
                    continue;

                node.isBound = false;

                m_nodesList.push_back(node);

                if(m_verboseOutput)
                {
                    std::cout << "Loading fluid node: " << "(";
                    for(unsigned short d = 0 ; d < m_dim ; ++d)
                    {
                        std::cout << node.position[d];
                        if(d == m_dim - 1)
                            std::cout << ")";
                        else
                            std::cout << ", ";
                    }
                    std::cout << std::endl;
                }
            }
        }
    }

#ifndef NDEBUG
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        for(std::size_t n2 = 0 ; n2 < m_nodesList.size() ; ++n2)
        {
            if(n != n2)
            {
                if(m_nodesList[n] == m_nodesList[n2])
                {
                    std::cerr << "Duplicate node found: " << "(";
                    for(unsigned short d = 0 ; d < m_dim ; ++d)
                    {
                        std::cerr << m_nodesList[n].position[d];
                        if(d == m_dim - 1)
                            std::cout << ")";
                        else
                            std::cout << ", ";
                    }
                    std::cerr << std::endl;
                    throw std::runtime_error("Duplicates nodes!");
                }
            }
        }
    }
#endif

    if(m_nodesList.empty())
        throw std::runtime_error("No nodes loaded! Did you add the physical groups in the .geo file?");

    gmsh::finalize();

    triangulateAlphaShape();
    //computeFreeSurfaceEdgeDetJ();
    computeElementsJ();
    computeElementsDetJ();
    computeElementsInvJ();
}

void Mesh::remesh()
{
    if(checkBoundingBox())
    {
        triangulateAlphaShape();
    }

    if(removeNodes())
    {
        triangulateAlphaShape();
        computeElementsJ();
        computeElementsDetJ();
    }

    addNodes();
    triangulateAlphaShape();
    //computeFreeSurfaceEdgeDetJ();
    computeElementsJ();
    computeElementsDetJ();
    computeElementsInvJ();
}

bool Mesh::removeNodes()
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh !");

    bool removeNodes = false;

    for(std::size_t i = 0 ; i < m_nodesList.size() ; ++i)
    {
        //We the node is free, it is not in an element; if a node is already tagged,
        //we can skip it (do not delete node twice)
        if(m_nodesList[i].toBeDeleted || m_nodesList[i].isFree)
            continue;

        for(unsigned int j = 0 ; j < m_nodesList[i].neighbourNodes.size() ; ++j)
        {
            double d = 0;

            for(unsigned short k = 0 ; k < m_dim ; ++k)
            {
                d +=(m_nodesList[i].position[k] - m_nodesList[m_nodesList[i].neighbourNodes[j]].position[k])
                   *(m_nodesList[i].position[k] - m_nodesList[m_nodesList[i].neighbourNodes[j]].position[k]);
            }
            d = std::sqrt(d);

            //Two nodes are too close.
            if(d <= m_gamma*m_hchar)
            {
                //If the neighbour nodes is touched, we delete the current nodes
                if(m_nodesList[m_nodesList[i].neighbourNodes[j]].touched)
                {
                    //Do not delete bounded or free surface nodes
                    if(m_nodesList[i].isBound || m_nodesList[i].isOnFreeSurface)
                        continue;

                    m_nodesList[i].toBeDeleted = true;
                    removeNodes = true;
                }
                //If the neighbour nodes is not touched, we delete the neoghbour nodes
                else
                {
                    //Do not delete bounded or free surface nodes
                    if(m_nodesList[m_nodesList[i].neighbourNodes[j]].isBound ||
                       m_nodesList[m_nodesList[i].neighbourNodes[j]].isOnFreeSurface)
                        continue;

                    m_nodesList[i].touched = true;
                    m_nodesList[m_nodesList[i].neighbourNodes[j]].toBeDeleted = true;
                    removeNodes = true;
                }
            }
        }
    }

    m_nodesList.erase(
    std::remove_if(m_nodesList.begin(), m_nodesList.end(), [this](const Node& node)
    {
       if(node.toBeDeleted)
       {
           if(this->m_verboseOutput)
           {
                std::cout << "Removing node " << "(";
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                {
                    std::cout << node.position[d];
                    if(d == m_dim - 1)
                        std::cout << ")";
                    else
                        std::cout << ", ";
                }
                std::cout << std::endl;
           }

           return true;
       }
       else
           return false;
    }), m_nodesList.end());

    return removeNodes;
}

void Mesh::restoreNodesList()
{
    if(m_nodesListSave.empty())
        throw std::runtime_error("The nodes list was not saved before or does not exist!");

    m_nodesList = m_nodesListSave;

    m_nodesListSave.clear();

    //computeFreeSurfaceEdgeDetJ();
    computeElementsJ();
    computeElementsDetJ();
    computeElementsInvJ();
}

void Mesh::saveNodesList()
{
    if(m_nodesList.empty())
        throw std::runtime_error("The nodes list does not exist!");

    m_nodesListSave = m_nodesList;
}

void Mesh::triangulateAlphaShape()
{
    if(m_dim == 2)
        triangulateAlphaShape2D();
    else if(m_dim == 3)
        triangulateAlphaShape3D();
    else
        throw std::runtime_error("Mesh dimension should be 2 or 3!");
}

void Mesh::updateNodesPosition(std::vector<double> deltaPos)
{
    if(deltaPos.size() != m_nodesList.size()*m_dim)
        throw std::runtime_error("Invalid size of the deltaPos vector");

    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        if(!m_nodesList[n].isBound)
        {
            for(unsigned short d = 0 ; d < m_dim ; ++d)
            {
                m_nodesList[n].position[d] += deltaPos[n + d*m_nodesList.size()];
            }
        }
    }

    //computeFreeSurfaceEdgeDetJ();
    computeElementsJ();
    computeElementsDetJ();
    computeElementsInvJ();
}

void Mesh::updateNodesPositionFromSave(std::vector<double> deltaPos)
{
    if(m_nodesListSave.empty())
        throw std::runtime_error("You did not save the nodes list!");
    else if(deltaPos.size() != m_nodesListSave.size()*m_dim)
        throw std::runtime_error("Invalid size of the deltaPos vector");

    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        if(!m_nodesList[n].isBound)
        {
            for(unsigned short d = 0 ; d < m_dim ; ++d)
            {
                m_nodesList[n].position[d] = m_nodesListSave[n].position[d] + deltaPos[n + d*m_nodesList.size()];
            }
        }
    }

    //computeFreeSurfaceEdgeDetJ();
    computeElementsJ();
    computeElementsDetJ();
    computeElementsInvJ();
}
