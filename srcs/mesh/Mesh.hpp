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

/**
 * \struct RemeshingParams
 * \brief Contains the parameters relative to the remeshing.
 */
struct RemeshingParams
{
    double hchar; /**< Characteristic size of an element (same as in .geo file). */
    double alpha; /**< Alpha parameter of the alpha-shape algorithm (triangles are
                       discared if  r_circumcircle > alpha*hchar). */
    double omega; /**< Control the addition of node if a triangle is too big (a node
                       is added if A_triangle > omege*hchar^2). */
    double gamma; /**< Control the deletetion of node if two are too close to each
                       other (a node is deleted if d_nodes < gamma*hchar). */
    std::vector<double> boundingBox; /**< Box delimiting the zone of nodes existence
                                          (format: [xmin, ymin, xmax, ymax]). */
};

/**
 * \class Mesh
 * \brief Represents a lagrangian mesh.
 */
class Mesh
{
    public:
        Mesh(const Params& params);
        ~Mesh();

        std::vector<Node> nodesList; /**< List of nodes of the mesh. */

        /**
         * \brief Add nodes in element whose area is too big (A_tringle > omega*hchar^2.
         * \return true if at least one node was added, false otherwise).
         */
        bool addNodes();

        /**
         * \brief Check if a node is outside the bounding box and deletes it if so.
         * \return true if at least one node was deleted, false otherwise.
         */
        bool checkBoundingBox();

        /**
         * \brief Compute the determinant of the of the Jacobian matrix for gauss
         *        integration for each triangle.
         */
        void computeDetJ();

        /**
         * \brief Compute the inverse Jacobian matrix for gauss integration for each
         *        triangle.
         */
        void computeInvJ();

        /**
         * \brief Get the determinant of the of the Jacobian matrix for gauss
         *        integration for a certain elm.
         * \param elm The element index.
         * \return The determinantf the of the Jacobian matrix
         */
        inline double getDetJ(std::size_t elm) const;

        /**
         * \brief Get the vector containing the index of the nodes in nodesList
                  making a certain element.
         * \param elm The index of the element.
         * \return The vector of index in nodesList.
         */
        inline std::vector<std::size_t> getElement(std::size_t elm) const;

        /**
         * \brief Return the number of element in the mesh.
         * \return The number of element in the mesh.
         */
        inline std::size_t getElementNumber() const;

        /**
         * \brief Get the ij element of the Jacobian matrix for a certain element.
         * \param elm The element index.
         * \param i The row index in the Jacobian matrix.
         * \param j The column index in the Jacobian matrix.
         * \return The ij element of the Jacobian matrix for a certain element.
         */
        inline double getInvJ(std::size_t elm, unsigned short i, unsigned short j) const;

        /**
         * \brief Get the mesh dimension (2 or 3).
         * \return The mesh dimension.
         */
        inline unsigned short getMeshDim() const;

        /**
         * \brief Get the gradient shape function matrix (linear shape functions,
         *        not multiplied by the Jacobian determinant) for a certain element.
         * \param elm The element index.
         * \return The gradient shape function matrix in the format:
         *         [dN1dx dN2dx dN3dx 0 0 0; 0 0 0 dN1dy dN2dy dN3dy; dN1dx dN2dx dN3dx dN1dy dN2dy dN3dy]
         */
        inline std::vector<Eigen::MatrixXd> getB(std::size_t elm) const;

        /**
         * \brief Get theshape function matrix (linear shape functions, not
         *        multiplied by the Jacobian determinant).
         * \return The gradient shape function matrix in the format:
         *         [N1 N2 N3 0 0 0; 0 0 0 N1 N2 N3]
         */
        inline std::vector<Eigen::MatrixXd> getN() const;

        /**
         * \brief Load the nodes from a file using gmsh.
         * \param fileName The name of the .msh file
         */
        void loadFromFile(std::string fileName);

        /**
         * \brief Remesh the nodes in nodesList using CGAL (Delaunay triangulation
         *        and alpha-shape)
         */
        void remesh();

        /**
         * \brief Removes nodes if they are too close from each other
         *       (d_nodes < gamma*hchar).
         * \return true if at least one node was deleted, false otherwise.
         */
        bool removeNodes();

    private:
         /**
         * \brief Compute the mesh dimension from the .msh file.
         * \return The mesh dimension (2 or 3).
         */
        inline unsigned short _computeMeshDim() const;

        RemeshingParams m_p; /**< Parameters of the remeshing. */

        bool m_verboseOutput; /**< Should the output be verbose? */

        std::vector<std::vector<std::size_t>> m_elementList; /**< The list of element (triplet of index in the nodesList. */
        std::vector<double> m_detJ; /**< The Jacobian matrix determinant of each element. */
        std::vector<Eigen::Matrix2d> m_invJ; /**< The inverse Jacobian matrix of each element. */

        unsigned short m_dim; /**< The mesh dimension. */

#ifdef DEBUG_GEOMVIEW
        CGAL::Geomview_stream m_gv; /**< A geomview stream object to interact with geomview if installed. */
#endif
};

#include "Mesh.inl"

#endif // MESH_HPP_INCLUDED
