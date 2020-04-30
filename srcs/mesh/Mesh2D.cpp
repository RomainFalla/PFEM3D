#include "Mesh.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
typedef Kernel::FT                                                          FT;
typedef CGAL::Triangulation_vertex_base_with_info_2<IndexType, Kernel>      Vb2;
typedef CGAL::Alpha_shape_vertex_base_2<Kernel, Vb2>                        asVb2;
typedef CGAL::Alpha_shape_face_base_2<Kernel>                               asFb2;
typedef CGAL::Triangulation_data_structure_2<asVb2,asFb2>                   asTds2;
typedef CGAL::Delaunay_triangulation_2<Kernel,asTds2>                       asTriangulation_2;
typedef CGAL::Alpha_shape_2<asTriangulation_2>                              Alpha_shape_2;
typedef Kernel::Point_2                                                     Point_2;
typedef Kernel::Triangle_2                                                  Triangle_2;


void Mesh::triangulateAlphaShape2D()
{
    if(m_nodesList.empty())
        throw std::runtime_error("You should load the mesh from a file before trying to remesh !");

    m_elementsList.clear();
//    m_freeSurfaceEdgesList.clear();

    // We have to construct an intermediate representation for CGAL. We also reset
    // nodes properties.
    std::vector<std::pair<Point_2, IndexType>> pointsList;
    for(IndexType i = 0 ; i < m_nodesList.size() ; ++i)
    {
        pointsList.push_back(std::make_pair(Point_2(m_nodesList[i].position[0],
                                                    m_nodesList[i].position[1]), i));

        m_nodesList[i].isFree = true;
        m_nodesList[i].isOnFreeSurface = false;
        m_nodesList[i].neighbourNodes.clear();
    }

    const Alpha_shape_2 as(pointsList.begin(), pointsList.end(),
                           FT(m_alpha*m_alpha*m_hchar*m_hchar),
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

            const std::vector<IndexType> element{face->vertex(0)->info(),
                                                 face->vertex(1)->info(),
                                                 face->vertex(2)->info()};

            // Those nodes are not free (flying nodes and not wetted boundary nodes)
            m_nodesList[element[0]].isFree = false;
            m_nodesList[element[1]].isFree = false;
            m_nodesList[element[2]].isFree = false;

            m_nodesList[element[0]].neighbourNodes.push_back(element[1]);
            m_nodesList[element[0]].neighbourNodes.push_back(element[2]);

            m_nodesList[element[1]].neighbourNodes.push_back(element[0]);
            m_nodesList[element[1]].neighbourNodes.push_back(element[2]);

            m_nodesList[element[2]].neighbourNodes.push_back(element[0]);
            m_nodesList[element[2]].neighbourNodes.push_back(element[1]);

            m_elementsList.push_back(element);
        }
    }

    for(auto& node : m_nodesList)
    {
        std::sort(node.neighbourNodes.begin(), node.neighbourNodes.end());
        node.neighbourNodes.erase(std::unique(node.neighbourNodes.begin(), node.neighbourNodes.end()), node.neighbourNodes.end());
    }

    for(Alpha_shape_2::Alpha_shape_edges_iterator it = as.alpha_shape_edges_begin() ;
        it != as.alpha_shape_edges_end() ; ++it)
    {
        // We compute the free surface nodes
        const Alpha_shape_2::Edge edgeAS {*it};
        if(as.classify(edgeAS) == Alpha_shape_2::REGULAR)
        {
            const std::vector<IndexType> edge{edgeAS.first->vertex((edgeAS.second+1)%3)->info(),
                                              edgeAS.first->vertex((edgeAS.second+2)%3)->info()};

            if(!(m_nodesList[edge[0]].isBound && m_nodesList[edge[1]].isBound))
            {
                m_nodesList[edge[0]].isOnFreeSurface = true;
                m_nodesList[edge[1]].isOnFreeSurface = true;

                //m_freeSurfaceEdgesList.push_back(edge);
            }
        }
    }

    // If an element is only composed of boundary nodes and the neighBournodes of
    // each of the three is equal to 2, this is a spurious triangle, we delete it
    m_elementsList.erase(
    std::remove_if(m_elementsList.begin(),  m_elementsList.end(), [this](const std::vector<IndexType>& element)
    {
        if(this->m_nodesList[element[0]].isBound && this->m_nodesList[element[1]].isBound &&
           this->m_nodesList[element[2]].isBound && this->m_nodesList[element[0]].neighbourNodes.size() == 2 &&
           this->m_nodesList[element[1]].neighbourNodes.size() == 2 && this->m_nodesList[element[2]].neighbourNodes.size() == 2)
        {
            this->m_nodesList[element[0]].isFree = true;
            this->m_nodesList[element[1]].isFree = true;
            this->m_nodesList[element[2]].isFree = true;

            return true;
        }
        else
            return false;

    }), m_elementsList.end());

    if(m_elementsList.empty())
        throw std::runtime_error("Something went wrong while remeshing. You might have not chosen a good \"hchar\" with regard to your .msh file");
}
