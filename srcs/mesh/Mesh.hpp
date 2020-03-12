#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include <vector>

#ifdef DEBUG_GEOMVIEW
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#endif
#include <gmsh.h>
#include <Eigen/Dense>

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
    double alpha; /**< Alpha parameter of the alpha-shape algorithm (triangles are discared if  r_circumcircle > alpha*hchar). */
    double omega; /**< Control the addition of node if a triangle is too big (a node is added if A_triangle > omege*hchar^2). */
    double gamma; /**< Control the deletetion of node if two are too close to each other (a node is deleted if d_nodes < gamma*hchar). */
    std::vector<double> boundingBox; /**< Box delimiting the zone of nodes existence (format: [xmin, ymin, xmax, ymax]). */
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

        /**
         * \brief Get the gradient shape function matrix (linear shape functions,
         *        not multiplied by the Jacobian determinant) for a certain element.
         * \param elm The element index.
         * \return The gradient shape function matrix in the format:
         *         [dN1dx dN2dx dN3dx 0 0 0; 0 0 0 dN1dy dN2dy dN3dy; dN1dx dN2dx dN3dx dN1dy dN2dy dN3dy]
         */
        inline std::vector<Eigen::MatrixXd> getB(std::size_t elm) const;

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
         * \brief Get theshape function matrix (linear shape functions, not
         *        multiplied by the Jacobian determinant).
         * \return The gradient shape function matrix in the format:
         *         [N1 N2 N3 0 0 0; 0 0 0 N1 N2 N3]
         */
        inline std::vector<Eigen::MatrixXd> getN() const;

        /**
         * \brief Get the number of nodes in the mesh.
         */
        inline std::size_t getNodesNumber() const;

        /**
         * \brief Get a coordinate of a node.
         * \param nodeIndex The index of the node in the nodes list.
         * \param coordinate The wanted coordinate (0, 1, ...).
         * \return The coordinate og the node.
         */
        inline double getNodePosition(std::size_t nodeIndex, unsigned short coordinate) const;

        /**
         * \brief Get a state of a node.
         * \param nodeIndex The index of the node in the nodes list.
         * \param state The wanted sate (0, 1, 2, 3, ...).
         * \return The value of the state for that node.
         */
        inline double getNodeState(std::size_t nodeIndex, unsigned short state) const;

        /**
         * \brief Check if a node is not attached to a fluid element.
         * \param nodeIndex The index of the node in the nodes list.
         * \return true is the node is free, false otherwize.
         */
        inline bool isNodeFree(std::size_t nodeIndex) const;

        /**
         * \brief Check if a node is a boundary node.
         * \param nodeIndex The index of the node in the nodes list.
         * \return true is the node is on the boundary, false otherwize.
         */
        inline bool isNodeBound(std::size_t nodeIndex) const;

        /**
         * \brief Check if a node is not attached to a fluid element.
         * \param nodeIndex The index of the node in the nodes list.
         * \return true is the node is a fluid input, false otherwize.
         */
        inline bool isNodeFluidInput(std::size_t nodeIndex) const;

        /**
         * \brief Load the nodes from a file using gmsh.
         * \param fileName The name of the .msh file
         */
        void loadFromFile(std::string fileName);

        /**
         * \brief Peform remeshing on the mesh.
         */
        void remesh();

        /**
         * \brief Restore the current nodes list from the list saved in saveNodesList.
         */
        void restoreNodesList();

        /**
         * \brief Save the current nodes list in another variable.
         */
        void saveNodesList();

        /**
         * \brief Update the nodes position.
         * \param deltaPos The variation of the coordinate;
         */
        void updateNodesPosition(std::vector<double> deltaPos);

        /**
         * \brief Update the nodes position (from the saved nodeS List).
         * \param deltaPos The variation of the coordinate;
         */
        void updateNodesPositionFromSave(std::vector<double> deltaPos);

        /**
         * \brief Set the position of a node.
         * \param nodeIndex The index of the node in the nodes list;
         * \param stateIndex The index of the state;
         * \param state The new value of the state;
         */
        inline void setNodeState(std::size_t nodeIndex, unsigned short stateIndex, double state);

        /**
         * \brief Set the number of states to be stored at node level.
         * \param statesNumber The number of state per nodes;
         */
        void setStatesNumber(unsigned short statesNumber);

    private:
        RemeshingParams m_p;    /**< Parameters of the remeshing. */
        bool m_verboseOutput;   /**< Should the output be verbose? */

        unsigned short m_dim;   /**< The mesh dimension. */

        std::vector<Node> m_nodesList;      /**< List of nodes of the mesh. */
        std::vector<Node> m_nodesListSave;  /**< A copy of the nodes list (usefull for non-linear algorithm). */
        std::vector<std::vector<std::size_t>> m_elementList;    /**< The list of element (triplet of index in the nodesList. */
        std::vector<double> m_detJ;         /**< The Jacobian matrix determinant of each element. */
        std::vector<std::vector<std::vector<double>>> m_invJ;   /**< The inverse Jacobian matrix of each element. */

#ifdef DEBUG_GEOMVIEW
        CGAL::Geomview_stream m_gv; /**< A geomview stream object to interact with geomview if installed. */
#endif

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
         * \brief Compute the mesh dimension from the .msh file.
         * \return The mesh dimension (2 or 3).
         */
        unsigned short computeMeshDim() const;

        /**
         * \brief Remesh the nodes in nodesList using CGAL (Delaunay triangulation
         *        and alpha-shape)
         */
        void triangulateAlphaShape();

        /**
         * \brief Removes nodes if they are too close from each other
         *       (d_nodes < gamma*hchar).
         * \return true if at least one node was deleted, false otherwise.
         */
        bool removeNodes();
};

#include "Mesh.inl"

#endif // MESH_HPP_INCLUDED
