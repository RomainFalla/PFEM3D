#include "Mesh.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h> 
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Regular_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
//typedef CGAL::Exact_predicates_exact_constructions_kernel                   Kernel;
typedef Kernel::FT                                                          FT;
typedef CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Kernel>    Vb2;
typedef CGAL::Alpha_shape_vertex_base_2<Kernel, Vb2>                        asVb2;
typedef CGAL::Alpha_shape_face_base_2<Kernel>                               asFb2;
typedef CGAL::Triangulation_data_structure_2<asVb2,asFb2>                   asTds2;
typedef CGAL::Delaunay_triangulation_2<Kernel,asTds2>                       asTriangulation_2;
typedef CGAL::Alpha_shape_2<asTriangulation_2>                              Alpha_shape_2;
typedef Kernel::Point_2                                                     Point_2;
typedef Kernel::Triangle_2                                                  Triangle_2;

// structure to keep track of the link between elements and nodes and transpose this part of cgal data structure into the code data structure.
struct FaceInfo2
{
    FaceInfo2() { isExt = true;  }
    std::size_t index;
    std::vector<std::size_t> neighboursElemIndexes;
    bool isExt;
};

typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,Kernel>        Fb2;
typedef CGAL::Triangulation_data_structure_2<Vb2, Fb2>                     Tds2;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds2>                       Triangulation_2;


//typedef CGAL::Regular_triangulation_vertex_base_2<Kernel, Vb2>              regVb2;
//typedef CGAL::Regular_triangulation_face_base_2<Kernel>                     regFb2;
//typedef CGAL::Alpha_shape_vertex_base_2<Kernel, regVb2>                     wasVb2;
//typedef CGAL::Alpha_shape_face_base_2<Kernel, regFb2>                       wasFb2;
//typedef CGAL::Triangulation_data_structure_2<wasVb2, wasFb2>                wasTds2;
//typedef CGAL::Regular_triangulation_2<Kernel, wasTds2>                      wasTriangulation_2;
//typedef CGAL::Alpha_shape_2<wasTriangulation_2>                             Weighted_Alpha_shape_2;
//typedef Kernel::Weighted_point_2                                            Weighted_point_2;



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
        if(as.classify(fit) == Alpha_shape_2::INTERIOR) // If true, the elements are fluid elements
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

            element.m_index = m_elementsList.size();

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

            facet.m_outNodeIndexes.push_back(edgeAS.first->vertex((edgeAS.second)%3)->info());

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

void Mesh::TriangulateWeightedAlphaShape2D() 
{
    if (m_nodesList.empty())
        throw std::runtime_error("You should load the mesh from a file before trying to remesh !");

    /*for (auto it = m_elementsList.begin(); it != m_elementsList.end(); ++it) 
    {
        Element el = *it;
        std::cout << "\n";
        std::cout << "Element:\n";
        for (auto it2 = el.m_nodesIndexes.begin(); it2 != el.m_nodesIndexes.end(); ++it2) 
        {
            std::cout << "w1 = " << m_nodesList[*it2].getWeight(m_alphaRatio, m_minTargetMeshSize) << "\n";
        }
    }*/
    //std::cout << "m_elementsList size = " << m_elementsList.size() << "\n";

    m_elementsList.clear();
    m_facetsList.clear();
    m_inBulkFacetList.clear();

    // We have to construct an intermediate representation for CGAL. We also reset
    // nodes properties.
    std::vector<std::pair<Point_2, std::size_t>> PointsList;

    for (std::size_t i = 0; i < m_nodesList.size(); ++i)
    {
        Point_2 P = Point_2(m_nodesList[i].m_position[0],m_nodesList[i].m_position[1]);
        //double weight = m_nodesList[i].getWeight(m_alphaRatio,m_minTargetMeshSize);
        //std::cout << "weight= " << weight << "\n";
        //Weighted_point_2 wP = Weighted_point_2(P, weight/100.);

        PointsList.push_back(std::make_pair(P, i));

        m_nodesList[i].m_neighbourNodes.clear();
        m_nodesList[i].m_elements.clear();
        m_nodesList[i].m_facets.clear();
    }

    //std::cout << "Weighted alpha shape ...!" << "\n";
    Triangulation_2 T(PointsList.begin(), PointsList.end());

    /*auto checkFaceDeletion = [&](Triangulation_2::Face_handle face) -> bool
    {
        std::size_t in0 = face->vertex(0)->info(), in1 = face->vertex(1)->info(), in2 = face->vertex(2)->info();

        if (m_nodesList[in0].isBound() && m_nodesList[in1].isBound() && m_nodesList[in2].isBound())
        {
            for (unsigned int i = 0; i <= 2; ++i)
            {
                std::set<Triangulation_2::Vertex_handle> neighbourVh;
                Triangulation_2::Face_circulator faceCirc = T.incident_faces(face->vertex(i)), done = faceCirc;
                do
                {
                    if (was.classify(faceCirc) == Weighted_Alpha_shape_2::INTERIOR)
                    {
                        for (unsigned int j = 0; j <= 2; ++j)
                        {
                            if (!was.is_infinite(faceCirc->vertex(j)))
                                neighbourVh.insert(faceCirc->vertex(j));
                        }
                    }
                    faceCirc++;
                } while (faceCirc != done);

                for (auto vh : neighbourVh)
                {
                    if (vh->info() == in0 || vh->info() == in1 || vh->info() == in2)
                        continue;

                    if (was.classify(vh) == Weighted_Alpha_shape_2::REGULAR || was.classify(vh) == Weighted_Alpha_shape_2::INTERIOR)
                    {
                        if (!m_nodesList[vh->info()].isBound())
                            return false;
                    }
                }
            }
            return true;
        }

        return false;
    };*/

    std::size_t index = 0;
    std::vector<Triangulation_2::Face_handle> elementCopies;
    std::vector<Triangulation_2::Face_handle> elToPotentiallyRemove;

    /*int count_total = 0;
    int count_bulk = 0;
    int count_inboundary = 0;
    int count_outboundary = 0;
    int count_outalpha = 0;
    int count_notSmooth = 0;*/

    //double meanElementSize = 0.;
    for (auto fit = T.finite_faces_begin(); fit != T.finite_faces_end(); ++fit)
    {
        //count_total++;
        index = m_elementsList.size();
        Triangulation_2::Face_handle face{ fit };
       
        std::size_t in0 = face->vertex(0)->info(), in1 = face->vertex(1)->info(), in2 = face->vertex(2)->info();
        std::vector<std::size_t> nodeIndexes = { in0,in1,in2 };
        ELEMENT_TYPE elmType = getElementType(nodeIndexes);

        Element element(*this);
        element.m_nodesIndexes = { in0, in1, in2 };

        element.build(nodeIndexes, m_nodesList);
        element.updateLargestExtension();

        double L = element.getLargestExtension();
        double minMeshSize = element.getMinMeshSize();
        double realMeshSize = element.getNaturalMeshSize(false);
        double naturalMeshSize = element.getNaturalMeshSize(true);
        double localMeshSize = element.getLocalMeshSize();

        double r = element.getCircumscribedRadius(); 
        //double r = getCircumScribedRadius(nodeIndexes);
        //meanElementSize += element.getSize();

        bool keepElement = false;

        double limitVal = 0.5 * m_minTargetMeshSize * m_minTargetMeshSize;

        if (L < 4. * m_maxProgressionFactor * naturalMeshSize)
        {
            std::vector<std::size_t> v = { in0,in1,in2 };
            
            if (elmType == IN_BULK)
                keepElement = true;
            else if (elmType == INSIDE_BOUNDARY)
            {
                //count_inboundary++;
                if (r > m_alphaRatio * localMeshSize)
                {
                    if (element.getSize() > limitVal)
                    {
                        keepElement = true;
                        double length = element.getBoundaryFacetSize();
                        if (realMeshSize/length < 0.5)
                            element.m_forcedRefinement = true;
                    }
                }
                else
                    keepElement = true;

            }
            else if (elmType == OUTSIDE_BOUNDARY && r < m_alphaRatio * localMeshSize && realMeshSize < 1.4 * minMeshSize)
                keepElement = true;

        }
        if (keepElement)
        {
            face->info().isExt = false;
            face->info().index = index;
            element.m_index = index;
            for (int j = 0; j < m_dim + 1; j++)
            {
                face->neighbor(j)->info().neighboursElemIndexes.push_back(index);
            }
            elementCopies.push_back(face);
            m_elementsList.push_back(std::move(element));
            for (std::size_t i : m_elementsList.back().m_nodesIndexes)
                m_nodesList[i].m_elements.push_back(m_elementsList.size() - 1);     
        }
    }
    std::size_t NbKeptElements = elementCopies.size();
    for (std::size_t i = 0; i < NbKeptElements; i++) 
    {
        m_elementsList[i].m_elemsIndexes = elementCopies[i]->info().neighboursElemIndexes; // copy the links from the CGAL data structure to our data structure.
    }

    elementCopies.clear();

#pragma omp parallel for default(shared)
    for (std::size_t n = 0; n < m_nodesList.size(); ++n)
    {
        std::sort(m_nodesList[n].m_neighbourNodes.begin(), m_nodesList[n].m_neighbourNodes.end());
        m_nodesList[n].m_neighbourNodes.erase(
            std::unique(m_nodesList[n].m_neighbourNodes.begin(), m_nodesList[n].m_neighbourNodes.end()),
            m_nodesList[n].m_neighbourNodes.end());
    }

    for (std::size_t i = 0; i < m_nodesList.size(); ++i)
    {
        if (m_nodesList[i].m_elements.size() != 0)
        {
            m_nodesList[i].m_isOnFreeSurface = false;
            m_nodesList[i].m_isOnBoundary = false;
        }
        else
        {   //flying nodes
            m_nodesList[i].m_isOnFreeSurface = true;
            m_nodesList[i].m_isOnBoundary = true;
        }
    }

    for (auto it = T.finite_edges_begin(); it != T.finite_edges_end(); ++it)
    {
        Triangulation_2::Edge edge1{ *it };
        Triangulation_2::Edge edge2=T.mirror_edge(edge1);

        Triangulation_2::Face_handle face1 = edge1.first;
        Triangulation_2::Face_handle face2 = edge2.first;

        Triangulation_2::Edge edgeToBuildFacet;

        Facet facet(*this);

        if (!face1->info().isExt && !face2->info().isExt) // not OnBoundary
        {
            facet.m_nodesIndexes = { edge1.first->vertex((edge1.second + 1) % 3)->info(),
                                    edge1.first->vertex((edge1.second + 2) % 3)->info() };

            facet.m_outNodeIndexes.push_back(edge1.first->vertex((edge1.second) % 3)->info());
            facet.m_outNodeIndexes.push_back(edge2.first->vertex((edge2.second) % 3)->info());
            facet.m_elementIndexes.push_back(edge1.first->info().index);
            facet.m_elementIndexes.push_back(edge2.first->info().index);
            m_inBulkFacetList.push_back(std::move(facet));
            continue;
        }

        if (face1->info().isExt && face2->info().isExt) // not on boundary + this condition --> not Regular
            continue;
       
        if (face1->info().isExt)
            edgeToBuildFacet = edge2;
        else 
            edgeToBuildFacet = edge1;

        facet.m_nodesIndexes = { edgeToBuildFacet.first->vertex((edgeToBuildFacet.second + 1) % 3)->info(),
                                edgeToBuildFacet.first->vertex((edgeToBuildFacet.second + 2) % 3)->info() };

        facet.m_outNodeIndexes.push_back(edgeToBuildFacet.first->vertex((edgeToBuildFacet.second) % 3)->info());
        facet.m_elementIndexes.push_back(edgeToBuildFacet.first->info().index);

        if (facet.isBound())
        {
            for (std::size_t k = 0; k < m_elementsList[facet.m_elementIndexes[0]].m_nodesIndexes.size(); k++)
            {
                if (facet.m_outNodeIndexes[0] == m_elementsList[facet.m_elementIndexes[0]].m_nodesIndexes[k])
                    m_elementsList[facet.m_elementIndexes[0]].m_chargedNode[k] = true;
            }
        }

        if (!(m_nodesList[facet.m_nodesIndexes[0]].isBound() && m_nodesList[facet.m_nodesIndexes[1]].isBound()))
        {
            m_nodesList[facet.m_nodesIndexes[0]].m_isOnFreeSurface = true;
            m_nodesList[facet.m_nodesIndexes[1]].m_isOnFreeSurface = true;
        }

        m_nodesList[facet.m_nodesIndexes[0]].m_isOnBoundary = true;
        m_nodesList[facet.m_nodesIndexes[1]].m_isOnBoundary = true;

        facet.computeJ();
        facet.computeDetJ();
        facet.computeInvJ();

        m_facetsList.push_back(std::move(facet));
    }

    std::size_t facetCounter = 0;
    for (auto it = m_inBulkFacetList.begin(); it != m_inBulkFacetList.end(); ++it)
    {
        m_elementsList[it->m_elementIndexes[0]].m_facets.push_back(facetCounter);
        m_elementsList[it->m_elementIndexes[1]].m_facets.push_back(facetCounter);
        for (std::size_t i : it->m_nodesIndexes)
            m_nodesList[i].m_facets.push_back(facetCounter); // because of this,  I have slightly modified the computation of the curvature function, cause we have to check the facets that are on the boundaries -R.F.
        facetCounter++;
    }
    for (auto it = m_facetsList.begin(); it != m_facetsList.end(); ++it)
    {
        m_elementsList[it->m_elementIndexes[0]].m_facets.push_back(facetCounter);
        for (std::size_t i : it->m_nodesIndexes)
            m_nodesList[i].m_facets.push_back(facetCounter);
        facetCounter++;
    }

    computeFSNormalCurvature();

    if (m_elementsList.empty())
        throw std::runtime_error("Something went wrong while remeshing. You might have not chosen a good \"hchar\" with regard to your .msh file");
        
}

void Mesh::computeFSNormalCurvature2D()
{
    if(!m_computeNormalCurvature)
        return;

    m_freeSurfaceCurvature.clear();
    m_boundFSNormal.clear();

    std::size_t minIndex = m_inBulkFacetList.size();
    for(std::size_t i = 0 ; i < m_nodesList.size() ; ++i)
    {
        const Node& node = m_nodesList[i];

        std::vector<Facet*> boundaryEdges;
        unsigned int nbFacets = node.getFacetCount();
     
        if(node.isOnFreeSurface() && !node.isBound())
        {
            assert(node.getFacetCount() != 0);

            for (unsigned int k = 0; k < nbFacets; k++) 
            {
                int index = int(node.m_facets[k] - minIndex);
                if (index >= 0)
                    boundaryEdges.push_back(&(m_facetsList[index]));
            }

            const Node& node_1 = ((boundaryEdges[0]->getNode(0) == node) ? boundaryEdges[0]->getNode(1) : boundaryEdges[0]->getNode(0));
            const Node& node1 = ((boundaryEdges[1]->getNode(0) == node) ? boundaryEdges[1]->getNode(1) : boundaryEdges[1]->getNode(0));

            //code unchanged from here...
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
            for (unsigned int k = 0; k < nbFacets; k++)
            {
                int index = int(node.m_facets[k] - minIndex);
                if (index >= 0)
                    boundaryEdges.push_back(&(m_facetsList[index]));
            }

            const Node& node_1 = ((boundaryEdges[0]->getNode(0) == node) ? boundaryEdges[0]->getNode(1) : boundaryEdges[0]->getNode(0));
            const Node& node1 = ((boundaryEdges[1]->getNode(0) == node) ? boundaryEdges[1]->getNode(1) : boundaryEdges[1]->getNode(0));

            const Node& outNode_1 = boundaryEdges[0]->getOutNode();
            const Node& outNode1 = boundaryEdges[1]->getOutNode();

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
