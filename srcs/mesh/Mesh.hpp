#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include <vector>

#ifdef DEBUG_GEOMVIEW
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#endif
#include <Eigen/Dense>
#include <gmsh.h>

#include "../quadrature/gausslegendre.hpp"
#include "Node.hpp"
#include "../params/Params.hpp"

struct RemeshingParams
{
    double hchar;
    double alpha;
    double omega;
    double gamma;
    std::vector<double> boundingBox;
};

class Mesh
{
    public:
        Mesh(const Params& params);
        ~Mesh();

        std::vector<Node> nodesList;

        bool addNodes();

        bool checkBoundingBox();

        void computeDetJ();
        void computeInvJ();

        inline double getDetJ(std::size_t elm) const;

        inline std::vector<std::size_t> getElement(std::size_t elm) const;
        inline std::size_t getElementNumber() const;

        inline double getInvJ(std::size_t elm, unsigned short i, unsigned short j) const;

        inline unsigned short getMeshDim() const;

        inline std::vector<Eigen::MatrixXd> getB(std::size_t elm) const;
        inline std::vector<Eigen::MatrixXd> getN() const;

        void loadFromFile(std::string fileName);

        void remesh();

        bool removeNodes();

    private:
        inline unsigned short _computeMeshDim() const;

        RemeshingParams m_p;

        bool m_verboseOutput;

        std::vector<std::vector<std::size_t>> m_elementList;
        std::vector<double> m_detJ;
        std::vector<Eigen::Matrix2d> m_invJ;

        unsigned short m_dim;

#ifdef DEBUG_GEOMVIEW
        CGAL::Geomview_stream m_gv;
#endif
};

#include "Mesh.inl"

#endif // MESH_HPP_INCLUDED
