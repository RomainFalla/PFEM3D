#include "Mesh.hpp"

#include <utility>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#include <Eigen/Dense>


typedef CGAL::Exact_predicates_inexact_constructions_kernel                     Kernel;
typedef Kernel::FT                                                              FT;
typedef CGAL::Triangulation_vertex_base_with_info_3<std::size_t, Kernel>        Vb3;
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
    m_facetsList.clear();

    // We have to construct an intermediate representation for CGAL. We also reset
    // nodes properties.
    std::vector<std::pair<Point_3, std::size_t>> pointsList;
    for(std::size_t i = 0 ; i < m_nodesList.size() ; ++i)
    {
        pointsList.push_back(std::make_pair(Point_3(m_nodesList[i].m_position[0],
                                                    m_nodesList[i].m_position[1],
                                                    m_nodesList[i].m_position[2]), i));

        m_nodesList[i].m_isOnFreeSurface = false;
        m_nodesList[i].m_neighbourNodes.clear();
        m_nodesList[i].m_elements.clear();
        m_nodesList[i].m_facets.clear();
    }

    const Alpha_shape_3 as(pointsList.begin(), pointsList.end(), m_alpha*m_alpha*m_hchar*m_hchar);

    auto checkCellDeletion = [&](Alpha_shape_3::Cell_handle cell) -> bool
    {
        std::size_t in0 = cell->vertex(0)->info(), in1 = cell->vertex(1)->info(),
                    in2 = cell->vertex(2)->info(), in3 = cell->vertex(3)->info();

        Tetrahedron_3 tetrahedron(pointsList[in0].first, pointsList[in1].first,
                                  pointsList[in2].first, pointsList[in3].first);

        if(tetrahedron.volume() < 1e-4*m_hchar*m_hchar*m_hchar)
            return true;

        if(m_nodesList[in0].isBound() && m_nodesList[in1].isBound() &&
           m_nodesList[in2].isBound() && m_nodesList[in3].isBound())
        {
            for(unsigned int i = 0 ; i < 4 ; ++i)
            {
                std::set<Alpha_shape_3::Vertex_handle> neighbourVh;
                as.adjacent_vertices(cell->vertex(i), std::inserter(neighbourVh, neighbourVh.begin()));

                for(auto vh: neighbourVh)
                {
                    if(as.is_infinite(vh))
                        continue;

                    if(vh->info() == in0 || vh->info() == in1 || vh->info() == in2 || vh->info() == in3)
                        continue;

                    if(as.classify(vh) == Alpha_shape_3::REGULAR || as.classify(vh) == Alpha_shape_3::INTERIOR)
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
    // perform operations on the remaining elements
    for(auto fit = as.finite_cells_begin() ; fit != as.finite_cells_end() ; ++fit)
    {
        // If true, the elements are fluid elements
        if(as.classify(fit) == Alpha_shape_3::INTERIOR)
        {
            const Alpha_shape_3::Cell_handle cell{fit};

            std::size_t in0 = cell->vertex(0)->info(), in1 = cell->vertex(1)->info(),
                        in2 = cell->vertex(2)->info(), in3 = cell->vertex(3)->info();

            if(checkCellDeletion(cell))
                continue;

            Element element(*this);
            element.m_nodesIndexes = {in0, in1, in2, in3};
            element.computeJ();
            element.computeDetJ();
            element.computeInvJ();

//            std::cout << in0 << ", " << in1 << ", " << in2 << ", " << in3 << ": " << m_nodesList.size() -1 << std::endl;
//            std::cout << m_nodesList[in0].m_neighbourNodes.size() << ", " << m_nodesList[in0].m_neighbourNodes.capacity() << std::endl;
            // We compute the neighbour nodes of each nodes
            m_nodesList[in0].m_neighbourNodes.push_back(in1);
            m_nodesList[in0].m_neighbourNodes.push_back(in2);
            m_nodesList[in0].m_neighbourNodes.push_back(in3);

            m_nodesList[in1].m_neighbourNodes.push_back(in0);
            m_nodesList[in1].m_neighbourNodes.push_back(in2);
            m_nodesList[in1].m_neighbourNodes.push_back(in3);

            m_nodesList[in2].m_neighbourNodes.push_back(in0);
            m_nodesList[in2].m_neighbourNodes.push_back(in1);
            m_nodesList[in2].m_neighbourNodes.push_back(in3);

            m_nodesList[in3].m_neighbourNodes.push_back(in0);
            m_nodesList[in3].m_neighbourNodes.push_back(in1);
            m_nodesList[in3].m_neighbourNodes.push_back(in2);

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

    for(auto it = as.finite_facets_begin() ; it != as.finite_facets_end() ; ++it)
    {
        // We compute the free surface nodes
        Alpha_shape_3::Facet facetAS {*it};
        if(as.classify(facetAS) == Alpha_shape_3::REGULAR)
        {
            Alpha_shape_3::Vertex_handle outVertex = facetAS.first->vertex((facetAS.second)%4);
            if(as.is_infinite(outVertex))
                facetAS = as.mirror_facet(facetAS);

            Alpha_shape_3::Cell_handle cell = facetAS.first;

            if(as.classify(cell) == Alpha_shape_3::EXTERIOR)
            {
                facetAS = as.mirror_facet(facetAS);
                cell = facetAS.first;
            }

            if(checkCellDeletion(cell))
                continue;

            Facet facet(*this);
            facet.m_nodesIndexes = {facetAS.first->vertex(asTriangulation_3::vertex_triple_index(facetAS.second, 0))->info(),
                                    facetAS.first->vertex(asTriangulation_3::vertex_triple_index(facetAS.second, 1))->info(),
                                    facetAS.first->vertex(asTriangulation_3::vertex_triple_index(facetAS.second, 2))->info()};
            facet.m_outNodeIndex = facetAS.first->vertex(facetAS.second)->info();
            facet.computeJ();
            facet.computeDetJ();
            facet.computeInvJ();

            if(facet.m_outNodeIndex == facet.m_nodesIndexes[0] ||
               facet.m_outNodeIndex == facet.m_nodesIndexes[1] ||
               facet.m_outNodeIndex == facet.m_nodesIndexes[2])
            {
                std::cerr << "A face with an outIndex equal to one of the face node index has been encountered" << std::endl;
                continue;
            }

            if(!(m_nodesList[facet.m_nodesIndexes[0]].isBound() &&
                 m_nodesList[facet.m_nodesIndexes[1]].isBound() &&
                 m_nodesList[facet.m_nodesIndexes[2]].isBound()))
            {
                m_nodesList[facet.m_nodesIndexes[0]].m_isOnFreeSurface = true;
                m_nodesList[facet.m_nodesIndexes[1]].m_isOnFreeSurface = true;
                m_nodesList[facet.m_nodesIndexes[2]].m_isOnFreeSurface = true;
            }

            m_facetsList.push_back(std::move(facet));

            for(std::size_t index: m_facetsList.back().m_nodesIndexes)
                m_nodesList[index].m_facets.push_back(m_facetsList.size() - 1);
        }
    }

    computeFSNormalCurvature3D();

    if(m_elementsList.empty())
        throw std::runtime_error("Something went wrong while remeshing. You might have not chosen a good \"hchar\" with regard to your .msh file");
}

void Mesh::regularTriangulateWeightedAlphaShape3D()
{
    if (m_nodesList.empty())
        throw std::runtime_error("You should load the mesh from a file before trying to remesh !");

    m_elementsList.clear();
    m_facetsList.clear();

    // We have to construct an intermediate representation for CGAL. We also reset
    // nodes properties.
    std::vector<std::pair<Point_3, std::size_t>> pointsList;
    for (std::size_t i = 0; i < m_nodesList.size(); ++i)
    {
        pointsList.push_back(std::make_pair(Point_3(m_nodesList[i].m_position[0],
            m_nodesList[i].m_position[1],
            m_nodesList[i].m_position[2]), i));

        m_nodesList[i].m_isOnFreeSurface = false;
        m_nodesList[i].m_neighbourNodes.clear();
        m_nodesList[i].m_elements.clear();
        m_nodesList[i].m_facets.clear();
    }

    const Alpha_shape_3 as(pointsList.begin(), pointsList.end(), m_alpha * m_alpha * m_hchar * m_hchar);

    auto checkCellDeletion = [&](Alpha_shape_3::Cell_handle cell) -> bool
    {
        std::size_t in0 = cell->vertex(0)->info(), in1 = cell->vertex(1)->info(),
            in2 = cell->vertex(2)->info(), in3 = cell->vertex(3)->info();

        Tetrahedron_3 tetrahedron(pointsList[in0].first, pointsList[in1].first,
            pointsList[in2].first, pointsList[in3].first);

        if (tetrahedron.volume() < 1e-4 * m_hchar * m_hchar * m_hchar)
            return true;

        if (m_nodesList[in0].isBound() && m_nodesList[in1].isBound() &&
            m_nodesList[in2].isBound() && m_nodesList[in3].isBound())
        {
            for (unsigned int i = 0; i < 4; ++i)
            {
                std::set<Alpha_shape_3::Vertex_handle> neighbourVh;
                as.adjacent_vertices(cell->vertex(i), std::inserter(neighbourVh, neighbourVh.begin()));

                for (auto vh : neighbourVh)
                {
                    if (as.is_infinite(vh))
                        continue;

                    if (vh->info() == in0 || vh->info() == in1 || vh->info() == in2 || vh->info() == in3)
                        continue;

                    if (as.classify(vh) == Alpha_shape_3::REGULAR || as.classify(vh) == Alpha_shape_3::INTERIOR)
                    {
                        if (!m_nodesList[vh->info()].isBound())
                            return false;
                    }
                }
            }

            return true;
        }

        return false;
    };

    // We check for each triangle which one will be kept (alpha shape), then we
    // perform operations on the remaining elements
    for (auto fit = as.finite_cells_begin(); fit != as.finite_cells_end(); ++fit)
    {
        // If true, the elements are fluid elements
        if (as.classify(fit) == Alpha_shape_3::INTERIOR)
        {
            const Alpha_shape_3::Cell_handle cell{ fit };

            std::size_t in0 = cell->vertex(0)->info(), in1 = cell->vertex(1)->info(),
                in2 = cell->vertex(2)->info(), in3 = cell->vertex(3)->info();

            if (checkCellDeletion(cell))
                continue;

            Element element(*this);
            element.m_nodesIndexes = { in0, in1, in2, in3 };
            element.computeJ();
            element.computeDetJ();
            element.computeInvJ();

            //            std::cout << in0 << ", " << in1 << ", " << in2 << ", " << in3 << ": " << m_nodesList.size() -1 << std::endl;
            //            std::cout << m_nodesList[in0].m_neighbourNodes.size() << ", " << m_nodesList[in0].m_neighbourNodes.capacity() << std::endl;
                        // We compute the neighbour nodes of each nodes
            m_nodesList[in0].m_neighbourNodes.push_back(in1);
            m_nodesList[in0].m_neighbourNodes.push_back(in2);
            m_nodesList[in0].m_neighbourNodes.push_back(in3);

            m_nodesList[in1].m_neighbourNodes.push_back(in0);
            m_nodesList[in1].m_neighbourNodes.push_back(in2);
            m_nodesList[in1].m_neighbourNodes.push_back(in3);

            m_nodesList[in2].m_neighbourNodes.push_back(in0);
            m_nodesList[in2].m_neighbourNodes.push_back(in1);
            m_nodesList[in2].m_neighbourNodes.push_back(in3);

            m_nodesList[in3].m_neighbourNodes.push_back(in0);
            m_nodesList[in3].m_neighbourNodes.push_back(in1);
            m_nodesList[in3].m_neighbourNodes.push_back(in2);

            m_elementsList.push_back(std::move(element));

            for (std::size_t index : m_elementsList.back().m_nodesIndexes)
                m_nodesList[index].m_elements.push_back(m_elementsList.size() - 1);
        }
    }

#pragma omp parallel for default(shared)
    for (std::size_t n = 0; n < m_nodesList.size(); ++n)
    {
        std::sort(m_nodesList[n].m_neighbourNodes.begin(), m_nodesList[n].m_neighbourNodes.end());
        m_nodesList[n].m_neighbourNodes.erase(
            std::unique(m_nodesList[n].m_neighbourNodes.begin(), m_nodesList[n].m_neighbourNodes.end()),
            m_nodesList[n].m_neighbourNodes.end());
    }

    for (auto it = as.finite_facets_begin(); it != as.finite_facets_end(); ++it)
    {
        // We compute the free surface nodes
        Alpha_shape_3::Facet facetAS{ *it };
        if (as.classify(facetAS) == Alpha_shape_3::REGULAR)
        {
            Alpha_shape_3::Vertex_handle outVertex = facetAS.first->vertex((facetAS.second) % 4);
            if (as.is_infinite(outVertex))
                facetAS = as.mirror_facet(facetAS);

            Alpha_shape_3::Cell_handle cell = facetAS.first;

            if (as.classify(cell) == Alpha_shape_3::EXTERIOR)
            {
                facetAS = as.mirror_facet(facetAS);
                cell = facetAS.first;
            }

            if (checkCellDeletion(cell))
                continue;

            Facet facet(*this);
            facet.m_nodesIndexes = { facetAS.first->vertex(asTriangulation_3::vertex_triple_index(facetAS.second, 0))->info(),
                                    facetAS.first->vertex(asTriangulation_3::vertex_triple_index(facetAS.second, 1))->info(),
                                    facetAS.first->vertex(asTriangulation_3::vertex_triple_index(facetAS.second, 2))->info() };
            facet.m_outNodeIndex = facetAS.first->vertex(facetAS.second)->info();
            facet.computeJ();
            facet.computeDetJ();
            facet.computeInvJ();

            if (facet.m_outNodeIndex == facet.m_nodesIndexes[0] ||
                facet.m_outNodeIndex == facet.m_nodesIndexes[1] ||
                facet.m_outNodeIndex == facet.m_nodesIndexes[2])
            {
                std::cerr << "A face with an outIndex equal to one of the face node index has been encountered" << std::endl;
                continue;
            }

            if (!(m_nodesList[facet.m_nodesIndexes[0]].isBound() &&
                m_nodesList[facet.m_nodesIndexes[1]].isBound() &&
                m_nodesList[facet.m_nodesIndexes[2]].isBound()))
            {
                m_nodesList[facet.m_nodesIndexes[0]].m_isOnFreeSurface = true;
                m_nodesList[facet.m_nodesIndexes[1]].m_isOnFreeSurface = true;
                m_nodesList[facet.m_nodesIndexes[2]].m_isOnFreeSurface = true;
            }

            m_facetsList.push_back(std::move(facet));

            for (std::size_t index : m_facetsList.back().m_nodesIndexes)
                m_nodesList[index].m_facets.push_back(m_facetsList.size() - 1);
        }
    }

    computeFSNormalCurvature3D();

    if (m_elementsList.empty())
        throw std::runtime_error("Something went wrong while remeshing. You might have not chosen a good \"hchar\" with regard to your .msh file");
}

void Mesh::computeFSNormalCurvature3D()
{
    if(!m_computeNormalCurvature)
        return;

    //TO DO: // this
    m_freeSurfaceCurvature.clear();
    m_boundFSNormal.clear();

    for(std::size_t i = 0 ; i < m_nodesList.size() ; ++i)
    {
        const Node& node = m_nodesList[i];
        if(!node.isOnFreeSurface() && !(node.isBound() && !node.isFree()))
            continue;

        assert(node.getFacetCount() != 0);

        //We first compute here a "to exterior" normal as the mean of faces normal
        std::array<double, 3> finalNodeNormal = {0, 0, 0};
        std::set<const Node*> pateletNodes;
        for(std::size_t j = 0; j < node.getFacetCount() ; ++j)
        {
            std::array<double, 3> facetNormal;

            const Facet& f = node.getFacet(j);

            const Node& outNode = f.getOutNode();

            std::array<const Node*, 2> facetNodes;

            if(f.getNode(0) == node)
            {
                facetNodes = {&f.getNode(1), &f.getNode(2)};
                pateletNodes.insert(&f.getNode(1));
                pateletNodes.insert(&f.getNode(2));
            }
            else if(f.getNode(1) == node)
            {
                facetNodes = {&f.getNode(0), &f.getNode(2)};
                pateletNodes.insert(&f.getNode(0));
                pateletNodes.insert(&f.getNode(2));
            }
            else
            {
                facetNodes = {&f.getNode(0), &f.getNode(1)};
                pateletNodes.insert(&f.getNode(0));
                pateletNodes.insert(&f.getNode(1));
            }

            const Node& n1 = *facetNodes[0];
            const Node& n2 = *facetNodes[1];

            std::array<double, 3> a = {
                n1.getCoordinate(0) - node.getCoordinate(0),
                n1.getCoordinate(1) - node.getCoordinate(1),
                n1.getCoordinate(2) - node.getCoordinate(2)
            };

            std::array<double, 3> b = {
                n2.getCoordinate(0) - node.getCoordinate(0),
                n2.getCoordinate(1) - node.getCoordinate(1),
                n2.getCoordinate(2) - node.getCoordinate(2)
            };

            facetNormal[0] = a[1]*b[2] - a[2]*b[1];
            facetNormal[1] = a[2]*b[0] - a[0]*b[2];
            facetNormal[2] = a[0]*b[1] - a[1]*b[0];

            double norm = std::sqrt(facetNormal[0]*facetNormal[0] + facetNormal[1]*facetNormal[1] + facetNormal[2]*facetNormal[2]);

            facetNormal[0] /= norm;
            facetNormal[1] /= norm;
            facetNormal[2] /= norm;

            double angle = std::acos((a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
                         /std::sqrt((a[0]*a[0]+a[1]*a[1]+a[2]*a[2])*(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])));

            std::array<double, 3> vecToOutNode = {
                outNode.getCoordinate(0) - node.getCoordinate(0),
                outNode.getCoordinate(1) - node.getCoordinate(1),
                outNode.getCoordinate(2) - node.getCoordinate(2)
            };

            if(vecToOutNode[0]*facetNormal[0] + vecToOutNode[1]*facetNormal[1] + vecToOutNode[2]*facetNormal[2] > 0)
            {
                facetNormal[0] *= -1.0;
                facetNormal[1] *= -1.0;
                facetNormal[2] *= -1.0;
            }

            for(std::size_t k = 0 ; k < finalNodeNormal.size() ; ++k)
                finalNodeNormal[k] += angle*facetNormal[k];
        }

        double finalNodeNormalNorm = std::sqrt(finalNodeNormal[0]*finalNodeNormal[0]
                                              + finalNodeNormal[1]*finalNodeNormal[1]
                                              + finalNodeNormal[2]*finalNodeNormal[2]);

        for(std::size_t k = 0 ; k < finalNodeNormal.size() ; ++k)
            finalNodeNormal[k] /= finalNodeNormalNorm;

        if(node.isBound())
            continue;

        //We now compute the two curvatures by fitting a 2D parabola onto the patelet of nodes.

        //Plane equation Ax+By+Cz+D = 0;
        double A = finalNodeNormal[0];
        double B = finalNodeNormal[1];
        double C = finalNodeNormal[2];
        double D = -(finalNodeNormal[0]*node.getCoordinate(0) + finalNodeNormal[1]*node.getCoordinate(1) + finalNodeNormal[2]*node.getCoordinate(2));

        // Determine two basis vector of the plane;
        std::array<double, 3> a, b;

        if(A != 0)
            a = {-(B + C)/A, 1, 1};

        else if(B != 0)
            a = {1, -(A + C)/B, 1};

        else if(C != 0)
            a = {1, 1, -(A + B)/C};

        else
            throw std::runtime_error("The computed normal was {0, 0, 0}");

        double aNorm = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

        for(std::size_t k = 0 ; k < finalNodeNormal.size() ; ++k)
            a[k] /= aNorm;

        b[0] = finalNodeNormal[1]*a[2] - finalNodeNormal[2]*a[1];
        b[1] = finalNodeNormal[2]*a[0] - finalNodeNormal[0]*a[2];
        b[2] = finalNodeNormal[0]*a[1] - finalNodeNormal[1]*a[0];

        //Determine in the local plane the coordinate of the patelet nodes
        std::vector<std::pair<double, double>> pateletLocalCoord;
        std::vector<double> pateletDist;
        for(const Node* pNode : pateletNodes)
        {
            double dist = A*pNode->getCoordinate(0) + B*pNode->getCoordinate(1) + C*pNode->getCoordinate(2) + D;

            pateletDist.push_back(dist);

            std::array<double, 3> distVec = {
                pNode->getCoordinate(0) - dist*finalNodeNormal[0] - node.getCoordinate(0),
                pNode->getCoordinate(1) - dist*finalNodeNormal[1] - node.getCoordinate(1),
                pNode->getCoordinate(2) - dist*finalNodeNormal[2] - node.getCoordinate(2)
            };

            pateletLocalCoord.push_back(std::make_pair(
                distVec[0]*a[0] + distVec[1]*a[1] + distVec[2]*a[2],
                distVec[0]*b[0] + distVec[1]*b[1] + distVec[2]*b[2]));
        }

        Eigen::VectorXd d = Eigen::VectorXd::Map(pateletDist.data(), pateletDist.size());

        Eigen::MatrixXd U(pateletDist.size(), 3);
        for(std::size_t k = 0 ; k < pateletNodes.size() ; ++k)
        {
            double u = pateletLocalCoord[k].first;
            double v = pateletLocalCoord[k].second;
            U(k, 0) = 0.5*u*u;
            U(k, 1) = u*v;
            U(k, 2) = 0.5*v*v;
        }

        Eigen::Matrix3d UTU = U.transpose() * U;
        Eigen::Vector3d UTd = U.transpose() * d;

        Eigen::Vector3d c = UTU.colPivHouseholderQr().solve(UTd);

        double delta = (c[0] - c[2])*(c[0] - c[2]) + 4*c[1]*c[1];
        double k1 = 0.5*(c[0] + c[2] + std::sqrt(delta));
        double k2 = 0.5*(c[0] + c[2] - std::sqrt(delta));

        m_freeSurfaceCurvature[i] = k1 + k2;
        m_boundFSNormal[i] = finalNodeNormal;
    }
}
