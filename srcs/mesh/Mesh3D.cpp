#include "Mesh.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel                     Kernel;
typedef Kernel::FT                                                              FT;
typedef CGAL::Triangulation_vertex_base_with_info_3<IndexType, Kernel>          Vb3;
typedef CGAL::Fixed_alpha_shape_vertex_base_3<Kernel, Vb3>                      asVb3;
typedef CGAL::Fixed_alpha_shape_cell_base_3<Kernel>                             asCb3;
typedef CGAL::Triangulation_data_structure_3<asVb3, asCb3>                      asTds3;
typedef CGAL::Delaunay_triangulation_3<Kernel, asTds3, CGAL::Fast_location>     asTriangulation_3;
typedef CGAL::Fixed_alpha_shape_3<asTriangulation_3>                            Alpha_shape_3;
typedef Kernel::Point_3                                                         Point_3;
typedef Kernel::Tetrahedron_3                                                   Tetrahedron_3;

void Mesh::triangulateAlphaShape3D()
{
    if(m_nodesList.empty())
        throw std::runtime_error("You should load the mesh from a file before trying to remesh !");

    m_elementsList.clear();
//    m_freeSurfaceEdgesList.clear();

    // We have to construct an intermediate representation for CGAL. We also reset
    // nodes properties.
    std::vector<std::pair<Point_3, IndexType>> pointsList;
    for(IndexType i = 0 ; i < m_nodesList.size() ; ++i)
    {
        pointsList.push_back(std::make_pair(Point_3(m_nodesList[i].position[0],
                                                    m_nodesList[i].position[1],
                                                    m_nodesList[i].position[2]), i));

        m_nodesList[i].isFree = true;
        m_nodesList[i].isOnFreeSurface = false;
        m_nodesList[i].neighbourNodes.clear();
        m_nodesList[i].belongingElements.clear();
    }

    const Alpha_shape_3 as(pointsList.begin(), pointsList.end(),
                           m_alpha*m_alpha*m_hchar*m_hchar);

    // We check for each triangle which one will be kept (alpha shape), then we
    // perfom operations on the remaining elements
    for(Alpha_shape_3::Finite_cells_iterator fit = as.finite_cells_begin() ;
        fit != as.finite_cells_end() ; ++fit)
    {
        // If true, the elements are fluid elements
        if(as.classify(fit) == Alpha_shape_3::INTERIOR)
        {
            const Alpha_shape_3::Cell_handle cell{fit};

            Element element = {};

            IndexType in0 = cell->vertex(0)->info();
            IndexType in1 = cell->vertex(1)->info();
            IndexType in2 = cell->vertex(2)->info();
            IndexType in3 = cell->vertex(3)->info();

            element.nodesIndexes = {in0, in1, in2, in3};

            Tetrahedron_3 tetrahedron(pointsList[in0].first,
                                      pointsList[in1].first,
                                      pointsList[in2].first,
                                      pointsList[in3].first);

            if(tetrahedron.volume() > 1e-4*m_hchar*m_hchar*m_hchar)
            {
                // Those nodes are not free (flying nodes and not wetted boundary nodes)
                m_nodesList[in0].isFree = false;
                m_nodesList[in1].isFree = false;
                m_nodesList[in2].isFree = false;
                m_nodesList[in3].isFree = false;

                // We compute the neighbour nodes of each nodes
                m_nodesList[in0].neighbourNodes.push_back(in1);
                m_nodesList[in0].neighbourNodes.push_back(in2);
                m_nodesList[in0].neighbourNodes.push_back(in3);

                m_nodesList[in1].neighbourNodes.push_back(in0);
                m_nodesList[in1].neighbourNodes.push_back(in2);
                m_nodesList[in1].neighbourNodes.push_back(in3);

                m_nodesList[in2].neighbourNodes.push_back(in0);
                m_nodesList[in2].neighbourNodes.push_back(in1);
                m_nodesList[in2].neighbourNodes.push_back(in3);

                m_nodesList[in3].neighbourNodes.push_back(in0);
                m_nodesList[in3].neighbourNodes.push_back(in1);
                m_nodesList[in3].neighbourNodes.push_back(in2);

                m_elementsList.push_back(std::move(element));

                m_nodesList[in0].belongingElements.push_back(m_elementsList.size() - 1);
                m_nodesList[in1].belongingElements.push_back(m_elementsList.size() - 1);
                m_nodesList[in2].belongingElements.push_back(m_elementsList.size() - 1);
                m_nodesList[in3].belongingElements.push_back(m_elementsList.size() - 1);
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(IndexType n = 0 ; n < m_nodesList.size() ; ++n)
    {
        std::sort(m_nodesList[n].neighbourNodes.begin(), m_nodesList[n].neighbourNodes.end());
        m_nodesList[n].neighbourNodes.erase(std::unique(m_nodesList[n].neighbourNodes.begin(), m_nodesList[n].neighbourNodes.end()), m_nodesList[n].neighbourNodes.end());
    }

    //This could take edges from bad elements-> construct tetrahedron from nodes edges to check
    for(Alpha_shape_3::Alpha_shape_facets_iterator it = as.alpha_shape_facets_begin() ;
        it != as.alpha_shape_facets_end() ; ++it)
    {
        // We compute the free surface nodes
        const Alpha_shape_3::Facet facetAS {*it};

        if(as.classify(facetAS) == Alpha_shape_3::REGULAR)
        {
            const std::vector<IndexType> edge{facetAS.first->vertex((facetAS.second+1)%4)->info(),
                                              facetAS.first->vertex((facetAS.second+2)%4)->info(),
                                              facetAS.first->vertex((facetAS.second+3)%4)->info()};

            if(m_nodesList[edge[0]].isBound && m_nodesList[edge[1]].isBound && m_nodesList[edge[2]].isBound)
                continue;

            m_nodesList[edge[0]].isOnFreeSurface = true;
            m_nodesList[edge[1]].isOnFreeSurface = true;
            m_nodesList[edge[2]].isOnFreeSurface = true;
                //m_freeSurfaceEdgesList.push_back(edge);
        }
    }

    // If an element is only composed of boundary nodes and the neighbour nodes of
    // each of the four are only boundary nodes, this is a spurious tetrahedron, we delete it
    m_elementsList.erase(
    std::remove_if(m_elementsList.begin(),  m_elementsList.end(),  [this](const Element& element)
    {
        if(this->m_nodesList[element.nodesIndexes[0]].isBound && this->m_nodesList[element.nodesIndexes[1]].isBound &&
           this->m_nodesList[element.nodesIndexes[2]].isBound && this->m_nodesList[element.nodesIndexes[3]].isBound)
        {
            for(unsigned short n = 0 ; n < element.nodesIndexes.size() ; ++n)
            {
                for(unsigned int i = 0 ; i < this->m_nodesList[element.nodesIndexes[n]].neighbourNodes.size() ; ++i)
                {
                    if(!this->m_nodesList[this->m_nodesList[element.nodesIndexes[n]].neighbourNodes[i]].isBound)
                    {
                        return false;
                    }
                }
            }

            this->m_nodesList[element.nodesIndexes[0]].isFree = true;
            this->m_nodesList[element.nodesIndexes[1]].isFree = true;
            this->m_nodesList[element.nodesIndexes[2]].isFree = true;
            this->m_nodesList[element.nodesIndexes[3]].isFree = true;

            return true;
        }
        else
            return false;

    }), m_elementsList.end());

    for(auto element : m_elementsList)
    {
        for(unsigned short n = 0 ; n < element.nodesIndexes.size() ; ++n)
            m_nodesList[element.nodesIndexes[n]].isFree = false;
    }

    #pragma omp parallel for default(shared)
    for(IndexType elm = 0 ; elm < m_elementsList.size() ; ++elm)
    {
        computeElementJ(m_elementsList[elm]);
        computeElementDetJ(m_elementsList[elm]);
        computeElementInvJ(m_elementsList[elm]);
    }

    if(m_elementsList.empty())
        throw std::runtime_error("Something went wrong while remeshing. You might have not chosen a good \"hchar\" with regard to your .msh file");
}
