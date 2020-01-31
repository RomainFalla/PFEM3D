#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include <vector>
#include <Eigen/Dense>
#include "../quadrature/gausslegendre.hpp"
#include "Node.hpp"
#include "../params/Params.hpp"

class Mesh
{
    public:
        Mesh(const Params& params);
        ~Mesh();

        std::vector<Node> nodesList;

        void addNodes();

        void computeDetInvJ();

        inline double getDetJ(std::size_t elm) const;

        inline std::vector<std::size_t> getElement(std::size_t elm) const;
        inline std::size_t getElementNumber() const;

        inline std::vector<bool> getIndices() const;

        inline double getInvJ(std::size_t elm, unsigned short i, unsigned short j) const;

        inline std::vector<Eigen::MatrixXd> getB(std::size_t elm) const;
        inline std::vector<Eigen::MatrixXd> getN() const;

        bool loadFromFile(std::string fileName);

        void remesh();

        void removeNodes();

    private:
        const Params& params;

        std::vector<std::vector<std::size_t>> elementList;
        std::vector<double> detJ;
        std::vector<Eigen::Matrix2d> invJ;
};

#include "Mesh.inl"

#endif // MESH_HPP_INCLUDED
