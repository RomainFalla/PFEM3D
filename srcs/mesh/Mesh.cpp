#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>

#ifdef DEBUG_GEOMVIEW
#include <chrono>
#include <thread>
#endif

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>

#include <nlohmann/json.hpp>

#include "Mesh.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Kernel> Vb;

typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds> Delaunay_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::FT FT;

typedef CGAL::Alpha_shape_vertex_base_2<Kernel, Vb>          asVb;
typedef CGAL::Alpha_shape_face_base_2<Kernel>                asFb;
typedef CGAL::Triangulation_data_structure_2<asVb,asFb>      asTds;
typedef CGAL::Delaunay_triangulation_2<Kernel,asTds>         asTriangulation_2;
typedef CGAL::Alpha_shape_2<asTriangulation_2>               Alpha_shape_2;


Mesh::Mesh(const Params& params)
#ifdef DEBUG_GEOMVIEW
: m_gv(CGAL::Bbox_3(0, 0, -10, 10, 10, 10))
#endif // DEBUG_GEOMVIEW
{
	nlohmann::json j = params.getJSON();

    m_verboseOutput = j["verboseOutput"].get<bool>();

	m_p.hchar = j["RemeshingParams"]["hchar"].get<double>();
	m_p.alpha = j["RemeshingParams"]["alpha"].get<double>();
	m_p.omega = j["RemeshingParams"]["omega"].get<double>();
	m_p.gamma = j["RemeshingParams"]["gamma"].get<double>();
	m_p.boundingBox = j["RemeshingParams"]["boundingBox"].get<std::vector<double>>();

	#ifdef DEBUG_GEOMVIEW
	m_gv.set_line_width(4);
    m_gv.set_bg_color(CGAL::Color(0, 200, 200));
	#endif // DEBUG_GEOMVIEW
}

Mesh::~Mesh()
{

}

bool Mesh::addNodes()
{
    assert(!m_elementList.empty() && !nodesList.empty() && "There is no mesh !");

    static const unsigned short dim = nodesList[0].position.size();
    static const unsigned short nUnknowns = nodesList[0].states.size();

    computeDetJ();

    bool addedNodes = false;

    for(std::size_t i = 0 ; i < m_elementList.size() ; ++i)
    {
        if(m_detJ[i]/2 > m_p.omega*m_p.hchar*m_p.hchar)
        {
            Node newNode(dim, nUnknowns);

            for(unsigned short k = 0 ; k < dim ; ++k)
            {
                newNode.position[k] = (nodesList[m_elementList[i][0]].position[k]
                                       + nodesList[m_elementList[i][1]].position[k]
                                       + nodesList[m_elementList[i][2]].position[k])/3.0;
            }

            for(unsigned short k = 0 ; k < nUnknowns ; ++k)
            {
                newNode.states[k] = (nodesList[m_elementList[i][0]].states[k]
                                     + nodesList[m_elementList[i][1]].states[k]
                                     + nodesList[m_elementList[i][2]].states[k])/3.0;
            }

            nodesList.push_back(newNode);

            if(m_verboseOutput)
            {
                std::cout << "Adding node (" << newNode.position[0] << ", "
                          << newNode.position[1] << ")" << std::endl;
            }

            addedNodes = true;
        }
    }

    return addedNodes;
}

bool Mesh::checkBoundingBox()
{
    assert(!m_elementList.empty() && !nodesList.empty() && "There is no mesh !");

    bool outofBBNodes = false;

    for(std::size_t n = 0 ; n < nodesList.size() ; ++n)
    {
        if(nodesList[n].position[0] < m_p.boundingBox[0]
           || nodesList[n].position[1] < m_p.boundingBox[1]
           || nodesList[n].position[0] > m_p.boundingBox[2]
           || nodesList[n].position[1] > m_p.boundingBox[3])
        {
            nodesList[n].toBeDeleted = true;
            outofBBNodes = true;
        }
    }

    nodesList.erase(std::remove_if(nodesList.begin(), nodesList.end(),
                                   [this](const Node& node){
                                       if(node.toBeDeleted)
                                       {
                                           if(this->m_verboseOutput)
                                           {
                                               std::cout << "Removing node out of bouding box ("
                                                         << node.position[0] << ", "
                                                         << node.position[1] << ") "
                                                         << std::endl;
                                           }

                                           return true;
                                       }
                                       else
                                           return false;
                                    }), nodesList.end());

    return outofBBNodes;
}

void Mesh::computeDetJ()
{
    assert(!m_elementList.empty() && !nodesList.empty() && "There is no mesh !");

    m_detJ.clear();
    m_detJ.resize(m_elementList.size());

    #pragma omp parallel for default(none) shared(m_detJ, m_elementList, nodesList)
    for(std::size_t i = 0 ; i < m_elementList.size() ; ++i)
    {
        m_detJ[i] = (nodesList[m_elementList[i][1]].position[0]
                         - nodesList[m_elementList[i][0]].position[0])
                         *(nodesList[m_elementList[i][2]].position[1]
                         - nodesList[m_elementList[i][0]].position[1])
                        - (nodesList[m_elementList[i][2]].position[0]
                         - nodesList[m_elementList[i][0]].position[0])
                         *(nodesList[m_elementList[i][1]].position[1]
                         - nodesList[m_elementList[i][0]].position[1]);
    }
}

void Mesh::computeInvJ()
{
    assert(!m_elementList.empty() && !nodesList.empty()
           && m_detJ.size() == m_elementList.size() && "There is no mesh or no detJ!");

    m_invJ.clear();
    m_invJ.resize(m_elementList.size());

    #pragma omp parallel for default(none) shared(m_invJ, m_detJ, m_elementList, nodesList)
    for(std::size_t i = 0 ; i < m_elementList.size() ; ++i)
    {
        m_invJ[i].resize(2,2);

        m_invJ[i](0,0) = nodesList[m_elementList[i][2]].position[1]
                             - nodesList[m_elementList[i][0]].position[1];

        m_invJ[i](0,1) = - nodesList[m_elementList[i][2]].position[0]
                             + nodesList[m_elementList[i][0]].position[0];

        m_invJ[i](1,0) = - nodesList[m_elementList[i][1]].position[1]
                             + nodesList[m_elementList[i][0]].position[1];

        m_invJ[i](1,1) = nodesList[m_elementList[i][1]].position[0]
                             - nodesList[m_elementList[i][0]].position[0];

        m_invJ[i] /= m_detJ[i];
    }
}

void Mesh::loadFromFile(std::string fileName)
{
    std::cout   << "================================================================"
                << std::endl
                << "                         LOADING THE MESH                       "
                << std::endl
                << "================================================================"
                << std::endl;

    nodesList.clear();

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);

    std::ifstream file(fileName);
    if(file.is_open())
        file.close();
    else
        throw std::runtime_error("The input .msh file does not exist!");

    gmsh::open(fileName);

    // Check that the mesh is not 3D
    m_dim = _computeMeshDim();
    if(m_dim != 2)
        throw std::runtime_error("Only 2D meshes supported currently!");

    if(m_p.boundingBox.size() != 2*m_dim)
        throw std::runtime_error("Bad bounding box format! Format: [xmin, ymin, xmax, ymax]");


    // We retrieve the tags of the physical groups of dimension m_dim and
    // m_dim-1
    std::vector<std::pair<int, int>> physGroupHandles2D;
    gmsh::model::getPhysicalGroups(physGroupHandles2D, m_dim);

    std::vector<std::pair<int, int>> physGroupHandles1D;
    gmsh::model::getPhysicalGroups(physGroupHandles1D, m_dim - 1);

    std::vector<std::size_t> nodesTagsBoundary;
    for(auto physGroup1D : physGroupHandles1D)
    {
        std::string name;
        gmsh::model::getPhysicalName(m_dim - 1, physGroup1D.second, name);
        if(name == "Boundary")
        {
            std::vector<double> coord;
            gmsh::model::mesh::getNodesForPhysicalGroup(m_dim - 1, physGroup1D.second,
                                                        nodesTagsBoundary, coord);

            for(std::size_t i = 0 ; i < nodesTagsBoundary.size() ; ++i)
            {
                Node node(m_dim, 3);
                node.position[0] = coord[3*i];
                node.position[1] = coord[3*i + 1];
                node.states[0] = 0;
                node.states[1] = 0;
                node.states[2] = 0;
                node.isBound = true;

                nodesList.push_back(node);

                if(m_verboseOutput)
                {
                    std::cout << "Loading boundary node: " << "(" << node.position[0]
                              << ", " << node.position[1] << ")" << std::endl;
                }
            }
        }
    }

    for(auto physGroup2D : physGroupHandles2D)
    {
        std::string name;
        gmsh::model::getPhysicalName(m_dim, physGroup2D.second, name);
        if(name == "Fluid")
        {
            std::vector<std::size_t> dummyNodesTags;
            std::vector<double> coord;
            gmsh::model::mesh::getNodesForPhysicalGroup(m_dim, physGroup2D.second,
                                                        dummyNodesTags, coord);

            for(std::size_t i = 0 ; i < dummyNodesTags.size() ; ++i)
            {
                //If the nodes is already on the boundary, we do not add it twice
                if(std::find(nodesTagsBoundary.begin(), nodesTagsBoundary.end(),
                             dummyNodesTags[i]) != nodesTagsBoundary.end())
                    continue;

                Node node(m_dim, 3);
                node.position[0] = coord[3*i];
                node.position[1] = coord[3*i + 1];
                node.states[0] = 0;
                node.states[1] = 0;
                node.states[2] = 0;
                node.isBound = false;

                nodesList.push_back(node);

                if(m_verboseOutput)
                {
                    std::cout << "Loading fluid node: " << "(" << node.position[0]
                              << ", " << node.position[1] << ")" << std::endl;
                }
            }
        }
    }

    gmsh::finalize();
}

void Mesh::remesh()
{
    if(nodesList.empty())
        throw std::runtime_error("You should load the mesh from a file before trying to remesh !");

    m_elementList.clear();

    // We have to construct an intermediate representation for CGAL. We also reset
    // nodes properties.
    std::vector<std::pair<Point_2, std::size_t>> pointsList;
    for(std::size_t i = 0 ; i < nodesList.size() ; ++i)
    {
        pointsList.push_back(std::make_pair(Point_2(nodesList[i].position[0],
                                                    nodesList[i].position[1]), i));

        nodesList[i].isFree = true;
        nodesList[i].isOnFreeSurface = false;
        nodesList[i].touched =  false;
        nodesList[i].toBeDeleted = false;
        nodesList[i].neighbourNodes.clear();
    }

    const Alpha_shape_2 as(pointsList.begin(), pointsList.end(),
                           FT(m_p.alpha*m_p.alpha*m_p.hchar*m_p.hchar),
                           Alpha_shape_2::GENERAL);

    // We check for each triangle which one will be kept (alpha shape), then we
    // perfom operations on the remaining elements
    for(Alpha_shape_2::Finite_faces_iterator fit = as.finite_faces_begin() ;
        fit != as.finite_faces_end() ; ++fit)
    {
        // If true, the elements are fluid elements
        if(as.classify(fit) == Alpha_shape_2::INTERIOR)
        {
            const Alpha_shape_2::Face_handle face{fit};

            const std::vector<std::size_t> element{face->vertex(0)->info(),
                                                   face->vertex(1)->info(),
                                                   face->vertex(2)->info()};

            // those nodes are not free (flying nodes and not wetted boundary nodes)
            nodesList[element[0]].isFree = false;
            nodesList[element[1]].isFree = false;
            nodesList[element[2]].isFree = false;

            // We compute the neighbour nodes of each nodes
            std::vector<std::size_t> temp;

            temp = nodesList[element[0]].neighbourNodes;
            if(std::find(temp.begin(), temp.end(), element[1]) == temp.end())
                nodesList[element[0]].neighbourNodes.push_back(element[1]);
            if(std::find(temp.begin(), temp.end(), element[2]) == temp.end())
                nodesList[element[0]].neighbourNodes.push_back(element[2]);

            temp = nodesList[element[1]].neighbourNodes;
            if(std::find(temp.begin(), temp.end(), element[0]) == temp.end())
                nodesList[element[1]].neighbourNodes.push_back(element[0]);
            if(std::find(temp.begin(), temp.end(), element[2]) == temp.end())
                nodesList[element[1]].neighbourNodes.push_back(element[2]);

            temp = nodesList[element[2]].neighbourNodes;
            if(std::find(temp.begin(), temp.end(), element[0]) == temp.end())
                nodesList[element[2]].neighbourNodes.push_back(element[0]);
            if(std::find(temp.begin(), temp.end(), element[1]) == temp.end())
                nodesList[element[2]].neighbourNodes.push_back(element[1]);

            m_elementList.push_back(element);
        }
    }

    for(Alpha_shape_2::Alpha_shape_vertices_iterator it = as.alpha_shape_vertices_begin() ;
        it != as.alpha_shape_vertices_end() ; ++it)
    {
        const Alpha_shape_2::Vertex_handle vert = *it;;
        if(!nodesList[vert->info()].isBound)
            nodesList[vert->info()].isOnFreeSurface = true;
    }

    // If an element is only composed of boundary nodes and the neighBournodes of
    // each of the three is equal to 2, this is a spurious triangle, we delete it
    m_elementList.erase(std::remove_if(m_elementList.begin(),  m_elementList.end(),
                            [this](const std::vector<std::size_t>& element)
                            {
                                if(this->nodesList[element[0]].isBound &&
                                   this->nodesList[element[1]].isBound &&
                                   this->nodesList[element[2]].isBound &&
                                   this->nodesList[element[0]].neighbourNodes.size() == 2 &&
                                   this->nodesList[element[1]].neighbourNodes.size() == 2 &&
                                   this->nodesList[element[2]].neighbourNodes.size() == 2)
                                {
                                    this->nodesList[element[0]].isFree = true;
                                    this->nodesList[element[1]].isFree = true;
                                    this->nodesList[element[2]].isFree = true;
                                    return true;
                                }
                                else
                                    return false;

                            }), m_elementList.end());

    if(m_elementList.empty())
        throw std::runtime_error("Something went wrong while remeshing. You might have not chosen a good \"hchar\" with regard to your .msh file");

#ifdef DEBUG_GEOMVIEW
    m_gv.clear();
    for(std::size_t n = 0 ; n < this->nodesList.size() ; ++n)
    {
        if(this->nodesList[n].isFree == false)
        {
            if(this->nodesList[n].isOnFreeSurface)
                m_gv << CGAL::Color(255, 0, 255);
            else
                m_gv << CGAL::Color(0, 0, 255);

            Point_2 point(this->nodesList[n].position[0],
                          this->nodesList[n].position[1]);

            m_gv << point;
        }
        else if(this->nodesList[n].isBound)
        {
            m_gv << CGAL::Color(255, 0, 0);
            Point_2 point(this->nodesList[n].position[0],
                          this->nodesList[n].position[1]);

            m_gv << point;
        }
        else
        {
            m_gv << CGAL::Color(0, 255, 0);
            Point_2 point(this->nodesList[n].position[0],
                          this->nodesList[n].position[1]);

            m_gv << point;
        }

    }

    m_gv.set_wired(true);
    m_gv << CGAL::Color(255, 255, 255);
    m_gv << as;

    std::this_thread::sleep_for(std::chrono::milliseconds(3000));
#endif // DEBUG_GEOMVIEW
}

bool Mesh::removeNodes()
{
    assert(!m_elementList.empty() && !nodesList.empty() && "There is no mesh !");

    bool removeNodes = false;

    for(std::size_t i = 0 ; i < nodesList.size() ; ++i)
    {
        if(nodesList[i].toBeDeleted || nodesList[i].isFree)
            continue;

        for(unsigned int j = 0 ; j < nodesList[i].neighbourNodes.size() ; ++j)
        {
            const double d = std::sqrt((nodesList[i].position[0]
                                     - nodesList[nodesList[i].neighbourNodes[j]].position[0])
                                     *(nodesList[i].position[0]
                                     - nodesList[nodesList[i].neighbourNodes[j]].position[0])
                                     +(nodesList[i].position[1]
                                     - nodesList[nodesList[i].neighbourNodes[j]].position[1])
                                     *(nodesList[i].position[1]
                                     - nodesList[nodesList[i].neighbourNodes[j]].position[1]));

            if(d <= m_p.gamma*m_p.hchar)
            {
                if(nodesList[nodesList[i].neighbourNodes[j]].touched)
                {
                    if(nodesList[i].isBound || nodesList[i].isOnFreeSurface)
                        continue;

                    nodesList[i].toBeDeleted = true;
                    removeNodes = true;
                }
                else
                {
                    if(nodesList[nodesList[i].neighbourNodes[j]].isBound ||
                       nodesList[nodesList[i].neighbourNodes[j]].isOnFreeSurface)
                        continue;

                    nodesList[i].touched = true;
                    nodesList[nodesList[i].neighbourNodes[j]].toBeDeleted = true;
                    removeNodes = true;
                }
            }
        }
    }

    nodesList.erase(std::remove_if(nodesList.begin(), nodesList.end(),
                                   [this](const Node& node){
                                       if(node.toBeDeleted)
                                       {
                                           if(this->m_verboseOutput)
                                           {
                                               std::cout << "Removing node ("
                                                         << node.position[0] << ", "
                                                         << node.position[1] << ") "
                                                         << std::endl;
                                           }

                                           return true;
                                       }
                                       else
                                           return false;
                                    }), nodesList.end());

    return removeNodes;
}
