#include "Mesh.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
typedef Kernel::FT                                                          FT;
typedef CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Kernel>    Vb2;
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
    m_facetsList.clear();

    // We have to construct an intermediate representation for CGAL. We also reset
    // nodes properties.
    std::vector<std::pair<Point_2, std::size_t>> pointsList;
    for(std::size_t i = 0 ; i < m_nodesList.size() ; ++i)
    {
        pointsList.push_back(std::make_pair(Point_2(m_nodesList[i].m_position[0],
                                                    m_nodesList[i].m_position[1]), i));

        m_nodesList[i].m_isOnFreeSurface = false;
        m_nodesList[i].m_neighbourNodes.clear();
        m_nodesList[i].m_elements.clear();
        m_nodesList[i].m_facets.clear();
    }

    const Alpha_shape_2 as(pointsList.begin(), pointsList.end(), m_alpha*m_alpha*m_hchar*m_hchar,
                           Alpha_shape_2::GENERAL);

    auto checkFaceDeletion = [&](Alpha_shape_2::Face_handle face) -> bool
    {
        std::size_t in0 = face->vertex(0)->info(), in1 = face->vertex(1)->info(), in2 = face->vertex(2)->info();

        if(m_nodesList[in0].isBound() && m_nodesList[in1].isBound() && m_nodesList[in2].isBound())
        {
            for(unsigned int i = 0 ; i <= 2 ; ++i)
            {
                std::set<Alpha_shape_2::Vertex_handle> neighbourVh;
                Alpha_shape_2::Face_circulator faceCirc = as.incident_faces(face->vertex(i)), done = faceCirc;
                do
                {
                    if(as.classify(faceCirc) == Alpha_shape_2::INTERIOR)
                    {
                        for(unsigned int j = 0 ; j <= 2 ; ++j)
                        {
                            if(!as.is_infinite(faceCirc->vertex(j)))
                                neighbourVh.insert(faceCirc->vertex(j));
                        }
                    }
                    faceCirc++;
                } while(faceCirc != done);

                for(auto vh : neighbourVh)
                {
                    if(vh->info() == in0 || vh->info() == in1 || vh->info() == in2)
                        continue;

                    if(as.classify(vh) == Alpha_shape_2::REGULAR || as.classify(vh) == Alpha_shape_2::INTERIOR)
                    {
                        if(!m_nodesList[vh->info()].isBound())
                            return false;
                    }
                }
            }
            return true;
        }

        return false;
    };

    // We check for each triangle which one will be kept (alpha shape), then we
    // perfom operations on the remaining elements
    for(auto fit = as.finite_faces_begin() ; fit != as.finite_faces_end() ; ++fit)
    {
        // If true, the elements are fluid elements
        if(as.classify(fit) == Alpha_shape_2::INTERIOR)
        {
            const Alpha_shape_2::Face_handle face{fit};
            std::size_t in0 = face->vertex(0)->info(), in1 = face->vertex(1)->info(), in2 = face->vertex(2)->info();

            if(checkFaceDeletion(face))
                continue;

            Element element(*this);
            element.m_nodesIndexes = {in0, in1, in2};

            element.computeJ();
            element.computeDetJ();
            element.computeInvJ();

            // Those nodes are not free (flying nodes and not wetted boundary nodes)
            m_nodesList[in0].m_neighbourNodes.push_back(in1);
            m_nodesList[in0].m_neighbourNodes.push_back(in2);

            m_nodesList[in1].m_neighbourNodes.push_back(in0);
            m_nodesList[in1].m_neighbourNodes.push_back(in2);

            m_nodesList[in2].m_neighbourNodes.push_back(in0);
            m_nodesList[in2].m_neighbourNodes.push_back(in1);

            m_elementsList.push_back(std::move(element));

            for(std::size_t index: m_elementsList.back().m_nodesIndexes)
                m_nodesList[index].m_elements.push_back(m_elementsList.size() - 1);
        }
    }

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        std::sort(m_nodesList[n].m_neighbourNodes.begin(), m_nodesList[n].m_neighbourNodes.end());
        m_nodesList[n].m_neighbourNodes.erase(
        std::unique(m_nodesList[n].m_neighbourNodes.begin(), m_nodesList[n].m_neighbourNodes.end()),
        m_nodesList[n].m_neighbourNodes.end());
    }

    for(auto it = as.finite_edges_begin() ; it != as.finite_edges_end() ; ++it)
    {
        // We compute the free surface nodes
        Alpha_shape_2::Edge edgeAS {*it};
        if(as.classify(edgeAS) == Alpha_shape_2::REGULAR)
        {
            Alpha_shape_2::Vertex_handle outVertex = edgeAS.first->vertex((edgeAS.second)%3);
            if(as.is_infinite(outVertex))
                edgeAS = as.mirror_edge(edgeAS);

            Alpha_shape_2::Face_handle face = edgeAS.first;

            if(as.classify(face) == Alpha_shape_2::EXTERIOR)
            {
                edgeAS = as.mirror_edge(edgeAS);
                face = edgeAS.first;
            }

            if(checkFaceDeletion(face))
                continue;

            Facet facet(*this);
            facet.m_nodesIndexes = {edgeAS.first->vertex((edgeAS.second+1)%3)->info(),
                                   edgeAS.first->vertex((edgeAS.second+2)%3)->info()};

            facet.m_outNodeIndex = edgeAS.first->vertex((edgeAS.second)%3)->info();

            if(!(m_nodesList[facet.m_nodesIndexes[0]].isBound() && m_nodesList[facet.m_nodesIndexes[1]].isBound()))
            {
                m_nodesList[facet.m_nodesIndexes[0]].m_isOnFreeSurface = true;
                m_nodesList[facet.m_nodesIndexes[1]].m_isOnFreeSurface = true;
            }

            facet.computeJ();
            facet.computeDetJ();
            facet.computeInvJ();

            m_facetsList.push_back(std::move(facet));

            for(std::size_t index: m_facetsList.back().m_nodesIndexes)
            m_nodesList[index].m_facets.push_back(m_facetsList.size() - 1);
        }
    }

    computeFSNormalCurvature();

    if(m_elementsList.empty())
        throw std::runtime_error("Something went wrong while remeshing. You might have not chosen a good \"hchar\" with regard to your .msh file");
}

void Mesh::computeFSNormalCurvature2D()
{
    if(!m_computeNormalCurvature)
        return;

    m_freeSurfaceCurvature.clear();
    m_boundFSNormal.clear();

    for(std::size_t i = 0 ; i < m_nodesList.size() ; ++i)
    {
        const Node& node = m_nodesList[i];
        if(node.isOnFreeSurface() && !node.isBound())
        {
            assert(node.getFacetCount() != 0);

            const Facet& f_1 = node.getFacet(0);
            const Facet& f1 = node.getFacet(1);

            const Node& node_1 = ((f_1.getNode(0) == node) ? f_1.getNode(1) : f_1.getNode(0));
            const Node& node1 = ((f1.getNode(0) == node) ? f1.getNode(1) : f1.getNode(0));

            double x_1 = node_1.getCoordinate(0);
            double x0 = node.getCoordinate(0);
            double x1 = node1.getCoordinate(0);

            double y_1 = node_1.getCoordinate(1);
            double y0 = node.getCoordinate(1);
            double y1 = node1.getCoordinate(1);

            double den = x0*(y1 - y_1) + x1*(y_1 - y0) + x_1*(y0 - y1);
            if(den == 0)
            {
                m_freeSurfaceCurvature[i] = 0;
                m_boundFSNormal[i] = {0, 0, 0};
                continue;
            }

            double m = - (x_1*x_1) - (y_1*y_1);
            double n = - (x0*x0) - (y0*y0);
            double o = - (x1*x1) - (y1*y1);

            double xc = -0.5*(m*(y0 - y1) + n*(y1 - y_1) - o*(y0-y_1))/den;
            double yc = -0.5*(-m*(x0 - x1) - n*(x1 - x_1) + o*(x0-x_1))/den;
            double rc = std::sqrt(xc*xc + yc*yc + (-m*(x0*y1 -x1*y0) -n*(x1*y_1 - x_1*y1) +o*(x0*y_1 - x_1*y0))/den);

            m_freeSurfaceCurvature[i] = 1/rc;

            std::array<double, 3> normal = {
                (xc - x0)/rc,
                (yc - y0)/rc,
                0
            };

            m_boundFSNormal[i] = normal;
        }
        else if(node.isBound() && !node.isFree())
        {
            const Facet& f_1 = node.getFacet(0);
            const Facet& f1 = node.getFacet(1);

            const Node& node_1 = ((f_1.getNode(0) == node) ? f_1.getNode(1) : f_1.getNode(0));
            const Node& node1 = ((f1.getNode(0) == node) ? f1.getNode(1) : f1.getNode(0));
            const Node& outNode_1 = f_1.getOutNode();
            const Node& outNode1 = f1.getOutNode();

            double x_1 = node_1.getCoordinate(0);
            double x0 = node.getCoordinate(0);
            double x1 = node1.getCoordinate(0);

            double y_1 = node_1.getCoordinate(1);
            double y0 = node.getCoordinate(1);
            double y1 = node1.getCoordinate(1);

            std::array<double, 3> normalF1 = {
                -(y0 - y1),
                x1 - x0,
                0
            };

            std::array<double, 3> vecToOutNode = {
                outNode1.getCoordinate(0) - x0,
                outNode1.getCoordinate(1) - y0,
                0
            };

            if(normalF1[0]*vecToOutNode[0] + normalF1[1]*vecToOutNode[1] > 0)
            {
                normalF1[0] *= -1;
                normalF1[1] *= -1;
            }

            double norm = std::sqrt(normalF1[0]*normalF1[0] + normalF1[1]*normalF1[1]);
            normalF1[0] /= norm;
            normalF1[1] /= norm;

            std::array<double, 3> normalF_1 = {
                -(y0 - y_1),
                x_1 - x0,
                0
            };

            vecToOutNode = {
                outNode_1.getCoordinate(0) - x0,
                outNode_1.getCoordinate(1) - y0,
                0
            };

            if(normalF_1[0]*vecToOutNode[0] + normalF_1[1]*vecToOutNode[1] > 0)
            {
                normalF_1[0] *= -1;
                normalF_1[1] *= -1;
            }

            norm = std::sqrt(normalF_1[0]*normalF_1[0] + normalF_1[1]*normalF_1[1]);
            normalF_1[0] /= norm;
            normalF_1[1] /= norm;

            std::array<double, 3> finalNodeNormal = {
                normalF1[0] + normalF_1[0],
                normalF1[1] + normalF_1[1],
                0
            };

            norm = std::sqrt(finalNodeNormal[0]*finalNodeNormal[0] + finalNodeNormal[1]*finalNodeNormal[1]);
            finalNodeNormal[0] /= norm;
            finalNodeNormal[1] /= norm;

            m_boundFSNormal[i] = finalNodeNormal;
        }
    }
}
