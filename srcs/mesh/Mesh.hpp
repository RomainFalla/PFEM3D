#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include <vector>
#include <Eigen/Dense>
#include "../quadrature/gausslegendre.hpp"
#include "Node.hpp"
#include "../params/Params.hpp"

#ifdef DEBUG_GEOMVIEW
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#endif

class Mesh
{
    public:
        Mesh(const Params& params);
        ~Mesh();

        std::vector<Node> nodesList;

        bool addNodes();

        void computeDetJ();
        void computeInvJ();

        inline double getDetJ(std::size_t elm) const;

        inline std::vector<std::size_t> getElement(std::size_t elm) const;
        inline std::size_t getElementNumber() const;

        inline double getInvJ(std::size_t elm, unsigned short i, unsigned short j) const;

        inline std::vector<Eigen::MatrixXd> getB(std::size_t elm) const;
        inline std::vector<Eigen::MatrixXd> getN() const;

        bool loadFromFile(std::string fileName);

        void remesh();

        bool removeNodes();

    private:
        const Params& m_params;

        std::vector<std::vector<std::size_t>> m_elementList;
        std::vector<double> m_detJ;
        std::vector<Eigen::Matrix2d> m_invJ;

#ifdef DEBUG_GEOMVIEW
        CGAL::Geomview_stream m_gv;
#endif
};

#include "Mesh.inl"

#endif // MESH_HPP_INCLUDED
