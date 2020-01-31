#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>

#ifdef DEBUG_GEOMVIEW
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#endif

#include "Mesh.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Kernel> Vb;

typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds> Delaunay_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::FT FT;

typedef CGAL::Alpha_shape_vertex_base_2<Kernel, Vb>          asVb;
typedef CGAL::Alpha_shape_face_base_2<Kernel>                asFb;
typedef CGAL::Triangulation_data_structure_2<asVb,asFb>      asTds;
typedef CGAL::Delaunay_triangulation_2<Kernel,asTds>         asTriangulation_2;
typedef CGAL::Alpha_shape_2<asTriangulation_2>               Alpha_shape_2;


Mesh::Mesh(const Params& params) :
params(params)
{

}

Mesh::~Mesh()
{

}

void Mesh::addNodes()
{
    assert(!this->elementList.empty() && !this->nodesList.empty() && "There is no mesh !");

    static const unsigned short dim = this->nodesList[0].position.size();
    static const unsigned short nUnknowns = this->nodesList[0].states.size();

    for(std::size_t i = 0 ; i < this->elementList.size() ; ++i)
    {
        if(this->detJ[i]/2 > this->params.omegaH2)
        {
            Node newNode(dim, nUnknowns);

            for(unsigned short k = 0 ; k < dim ; ++k)
            {
                newNode.position[k] = (1/3)*(this->nodesList[this->elementList[i][0]].position[k]
                                             + this->nodesList[this->elementList[i][1]].position[k]
                                             + this->nodesList[this->elementList[i][2]].position[k]);
            }

            for(unsigned short k = 0 ; k < nUnknowns ; ++k)
            {
                newNode.states[k] = (1/3)*(this->nodesList[this->elementList[i][0]].states[k]
                                           + this->nodesList[this->elementList[i][1]].states[k]
                                           + this->nodesList[this->elementList[i][2]].states[k]);
            }

            this->nodesList.push_back(newNode);
        }
    }
}

void Mesh::computeDetInvJ()
{
    assert(!this->elementList.empty() && !this->nodesList.empty() && "There is no mesh !");

    this->detJ.clear();
    this->invJ.clear();
    this->detJ.resize(this->elementList.size());
    this->invJ.resize(this->elementList.size());

    for(std::size_t i = 0 ; i < this->elementList.size() ; ++i)
    {
        this->detJ[i] = (this->nodesList[this->elementList[i][1]].position[0]
                         - this->nodesList[this->elementList[i][0]].position[0])
                         *(this->nodesList[this->elementList[i][2]].position[1]
                         - this->nodesList[this->elementList[i][0]].position[1])
                        - (this->nodesList[this->elementList[i][2]].position[0]
                         - this->nodesList[this->elementList[i][0]].position[0])
                         *(this->nodesList[this->elementList[i][1]].position[1]
                         - this->nodesList[this->elementList[i][0]].position[1]);

        this->invJ[i].resize(2,2);

        this->invJ[i](0,0) = this->nodesList[this->elementList[i][2]].position[1]
                             - this->nodesList[this->elementList[i][0]].position[1];

        this->invJ[i](0,1) = - this->nodesList[this->elementList[i][2]].position[0]
                             + this->nodesList[this->elementList[i][0]].position[0];

        this->invJ[i](1,0) = - this->nodesList[this->elementList[i][1]].position[1]
                             + this->nodesList[this->elementList[i][0]].position[1];

        this->invJ[i](1,1) = this->nodesList[this->elementList[i][1]].position[0]
                             - this->nodesList[this->elementList[i][0]].position[0];

        this->invJ[i] /= this->detJ[i];
    }
}

bool Mesh::loadFromFile(std::string fileName)
{
    this->nodesList.clear();

    for (unsigned int i = 1 ; i <= static_cast<unsigned int>(5.0/params.hchar) ; ++i)
    {
        for (unsigned int j = 1 ; j <= static_cast<unsigned int>(5.0/params.hchar) ; ++j)
        {
            Node node(2, 3);
            node.position[0] = i*params.hchar;
            node.position[1] = j*params.hchar;
            node.states[0] = 0;
            node.states[1] = 0;
            node.states[2] = (5-node.position[i])*params.gravity*params.fluidParameters[0];
            node.isBound = false;
            this->nodesList.push_back(node);
        }
    }

    for (unsigned int i = 0 ; i < static_cast<unsigned int>(10.0/params.hchar) ; ++i)
    {
        Node node(2, 3);
        node.position[0] = 0*params.hchar;
        node.position[1] = i*params.hchar;
        node.states[0] = 0;
        node.states[1] = 0;
        node.states[2] = 0*(5-node.position[i])*params.gravity*params.fluidParameters[0];
        node.isBound = true;
        nodesList.push_back(node);
        node.position[0] = 10;
        node.position[1] = i*params.hchar;
        node.isBound = true;
        this->nodesList.push_back(node);
        if(i != 0)
        {
            node.position[0] = i*params.hchar;
            node.position[1] = 0;
            node.states[0] = 0;
            node.states[1] = 0;
            node.states[2] = 0*(5-node.position[i])*params.gravity*params.fluidParameters[0];
            node.isBound = true;
            this->nodesList.push_back(node);
        }
    }

    return true;
}

void Mesh::remesh()
{
#ifdef DEBUG_GEOMVIEW
    CGAL::Geomview_stream gv(CGAL::Bbox_3(0, 0, -10, 10, 10, 10));
    gv.set_line_width(4);
    gv.set_bg_color(CGAL::Color(0, 200, 200));
#endif

    if(this->nodesList.empty())
        throw std::runtime_error("You should load the mesh from a file before trying to remesh !");

    this->elementList.clear();

    // We have to construct an intermediate representation for CGAL. We also reset
    // nodes properties.
    std::vector<std::pair<Point_2, std::size_t>> pointsList;
    for(std::size_t i = 0 ; i < this->nodesList.size() ; ++i)
    {
        pointsList.push_back(std::make_pair(Point_2(this->nodesList[i].position[0],
                                                    this->nodesList[i].position[1]), i));

        this->nodesList[i].isFree = true;
        this->nodesList[i].isOnFreeSurface = false;
        this->nodesList[i].touched=  false;
        this->nodesList[i].toBeDeleted = false;
        this->nodesList[i].neighbourNodes.clear();
    }

    const Alpha_shape_2 as(pointsList.begin(), pointsList.end(),
                           FT(this->params.alphaHchar*this->params.alphaHchar),
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
#ifdef DEBUG_GEOMVIEW
            gv << CGAL::Color(0, 0, 255);
            gv << face->vertex(0)->point();
            gv << face->vertex(1)->point();
            gv << face->vertex(2)->point();
#endif

            this->elementList.push_back(element);

            // those nodes are not free (flying nodes and not wetted boundary nodes)
            this->nodesList[element[0]].isFree = false;
            this->nodesList[element[1]].isFree = false;
            this->nodesList[element[2]].isFree = false;

            // We compute the neighbour nodes of each nodes
            std::vector<std::size_t> temp;

            temp = this->nodesList[element[0]].neighbourNodes;
            if(std::find(temp.begin(), temp.end(), element[1]) == temp.end())
                this->nodesList[element[0]].neighbourNodes.push_back(element[1]);
            if(std::find(temp.begin(), temp.end(), element[2]) == temp.end())
                this->nodesList[element[0]].neighbourNodes.push_back(element[2]);

            temp = this->nodesList[element[1]].neighbourNodes;
            if(std::find(temp.begin(), temp.end(), element[0]) == temp.end())
                this->nodesList[element[1]].neighbourNodes.push_back(element[0]);
            if(std::find(temp.begin(), temp.end(), element[2]) == temp.end())
                this->nodesList[element[1]].neighbourNodes.push_back(element[2]);

            temp = nodesList[element[2]].neighbourNodes;
            if(std::find(temp.begin(), temp.end(), element[0]) == temp.end())
                this->nodesList[element[2]].neighbourNodes.push_back(element[0]);
            if(std::find(temp.begin(), temp.end(), element[1]) == temp.end())
                this->nodesList[element[2]].neighbourNodes.push_back(element[1]);
        }
    }

    for(Alpha_shape_2::Alpha_shape_vertices_iterator it = as.alpha_shape_vertices_begin() ;
        it != as.alpha_shape_vertices_end() ; ++it)
    {
        const Alpha_shape_2::Vertex_handle vert = *it;;
        if(!this->nodesList[vert->info()].isBound)
            this->nodesList[vert->info()].isOnFreeSurface = true;
    }

    // If an element is only composed of boundary nodes and the neighBournodes of
    // each of the three is equal to 2, this is a spurious triangle, we delete it
    elementList.erase(std::remove_if(elementList.begin(),  elementList.end(),
                            [this](const std::vector<std::size_t>& element)
                            {
                                if(this->nodesList[element[0]].isBound &&
                                   this->nodesList[element[1]].isBound &&
                                   this->nodesList[element[2]].isBound &&
                                   this->nodesList[element[0]].neighbourNodes.size() == 2 &&
                                   this->nodesList[element[1]].neighbourNodes.size() == 2 &&
                                   this->nodesList[element[2]].neighbourNodes.size() == 2)
                                    return true;
                                else
                                    return false;

                            }), elementList.end());
#ifdef DEBUG_GEOMVIEW
    gv.set_wired(true);
    gv << CGAL::Color(255, 255, 255);
    gv << as;

    std::cout << "Enter a key to finish" << std::endl;
    char ch2;
    std::cin >> ch2;
#endif
}

void Mesh::removeNodes()
{
    assert(!this->elementList.empty() && !this->nodesList.empty() && "There is no mesh !");

    for(std::size_t i = 0 ; i < this->nodesList.size() ; ++i)
    {
        if(this->nodesList[i].toBeDeleted || this->nodesList[i].isFree)
            continue;

        for(unsigned int j = 0 ; this->nodesList[i].neighbourNodes.size() ; ++j)
        {
            const double d = std::sqrt((this->nodesList[i].position[0]
                                     - this->nodesList[this->nodesList[i].neighbourNodes[j]].position[0])
                                     *(this->nodesList[i].position[0]
                                     - this->nodesList[this->nodesList[i].neighbourNodes[j]].position[0])
                                     +(this->nodesList[i].position[1]
                                     - this->nodesList[this->nodesList[i].neighbourNodes[j]].position[1])
                                     *(this->nodesList[i].position[1]
                                     - this->nodesList[this->nodesList[i].neighbourNodes[j]].position[1]));

            if(d <= this->params.gammaH)
            {
                if(this->nodesList[this->nodesList[i].neighbourNodes[j]].touched)
                {
                    if(this->nodesList[i].isBound || this->nodesList[i].isOnFreeSurface)
                        continue;

                    this->nodesList[i].toBeDeleted = true;
                }
                else
                {
                    if(this->nodesList[this->nodesList[i].neighbourNodes[j]].isBound ||
                       this->nodesList[this->nodesList[i].neighbourNodes[j]].isOnFreeSurface)
                        continue;

                    this->nodesList[i].touched = true;
                    this->nodesList[this->nodesList[i].neighbourNodes[j]].toBeDeleted = true;
                }
            }
        }
    }

    this->nodesList.erase(std::remove_if(this->nodesList.begin(), this->nodesList.end(),
                                         [](const Node& node){
                                            if(node.toBeDeleted)
                                                return true;
                                            else
                                                return false;
                                         }), this->nodesList.end());
}
