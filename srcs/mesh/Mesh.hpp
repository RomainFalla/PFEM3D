#pragma once
#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include <nlohmann/json.hpp>

#include "Node.hpp"

#include "Mesh_export.h"

/**
 * \class Mesh
 * \brief Represents a lagrangian mesh.
 */
class MESH_API Mesh
{
    public:
        Mesh(const nlohmann::json& j);
        ~Mesh();

        /**
         * \return The characteristic alpha-shape parameter.
         */
        inline double getAlpha() const;

        /**
         * \param elm The index of the element.
         * \return The vector containing the index of the nodes in nodesList
         *         making a certain element.
         */
        inline std::vector<std::size_t> getElement(std::size_t elm) const;

        /**
         * \param elm The element index.
         * \return The determinant of the of the Jacobian matrix
         */
        inline double getElementDetJ(std::size_t elm) const;

        /**
         * \param elm The element index.
         * \param i The row index in the Jacobian matrix.
         * \param j The column index in the Jacobian matrix.
         * \return The ij element of the Jacobian matrix for a certain element.
         */
        inline double getElementInvJ(std::size_t elm, unsigned short i, unsigned short j) const;

        /**
         * \return The number of elements in the mesh.
         */
        inline std::size_t getElementsNumber() const;

//        /**
//         * \param elm The element index.
//         * \return The determinant of the of the Jacobian matrix
//         */
//        inline double getFreeSurfaceDetJ(std::size_t edge) const;
//
//
//        /**
//         * \param edge The index of the edge.
//         * \return The vector containing the index of the nodes in nodesList
//         *         making a certain edge.
//         */
//        inline std::vector<std::size_t> getFreeSurfaceEdge(std::size_t edge) const;
//
//        /**
//         * \return The number of free surface edges in the mesh.
//         */
//        inline std::size_t getFreeSurfaceEdgesNumber() const;

        /**
         * \return The characteristic node removing parameter.
         */
        inline double getGamma() const;

        /**
         * \param point The index to the considered Gauss point.
         * \param coordinate The index of the coordinate we want.
         * \return The coordinate of the wanted Gauss point.
         */
        inline double getGaussPoints(unsigned short point, unsigned short coordinate) const;

        /**
         * \return The number of gauss point used (depends on the mesh dimension).
         */
        inline unsigned short getGaussPointsNumber() const;

        /**
         * \param point The index to the considered Gauss point.
         * \return The weight associated to this Gauss Point.
         */
        inline double getGaussWeight(unsigned short point) const;

        /**
         * \return The characteristic size of the mesh.
         */
        inline double getHchar() const;

        /**
         * \return The mesh dimension.
         */
        inline unsigned short getMeshDim() const;

        /**
         * \return The number of nodes in the mesh.
         */
        inline std::size_t getNodesNumber() const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \param coordinate The wanted coordinate (0, 1, ...).
         * \return The coordinate of the node.
         */
        inline double getNodePosition(std::size_t nodeIndex, unsigned short coordinate) const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \param state The wanted sate (0, 1, 2, 3, ...).
         * \return The value of the state for that node.
         */
        inline double getNodeState(std::size_t nodeIndex, unsigned short state) const;

        /**
         * \return The characteristic node adding parameter.
         */
        inline double getOmega() const;

        /**
         * \return The size of the element in the reference coordinate system.
         */
        inline double getRefElementSize() const;

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
        void loadFromFile(const std::string& fileName);

        /**
         * \brief Perform remeshing on the mesh.
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
        inline void setStatesNumber(unsigned short statesNumber);

    private:
        bool m_verboseOutput;   /**< Should the output be verbose? */

        double m_hchar; /**< Characteristic size of an element (same as in .geo file). */
        double m_alpha; /**< Alpha parameter of the alpha-shape algorithm (triangles are discared if  r_circumcircle > alpha*hchar). */
        double m_omega; /**< Control the addition of node if a triangle is too big (a node is added if A_triangle > omege*hchar^2). */
        double m_gamma; /**< Control the deletetion of node if two are too close to each other (a node is deleted if d_nodes < gamma*hchar). */
        std::vector<double> m_boundingBox; /**< Box delimiting the zone of nodes existence (format: [xmin, ymin, xmax, ymax]). */

        unsigned short m_dim;   /**< The mesh dimension. */

        std::vector<Node> m_nodesList;      /**< List of nodes of the mesh. */
        std::vector<Node> m_nodesListSave;  /**< A copy of the nodes list (usefull for non-linear algorithm). */
        std::vector<std::vector<std::size_t>> m_elementsList;    /**< The list of element (triplet of index in the nodesList. */
        //std::vector<std::vector<std::size_t>> m_freeSurfaceEdgesList;   /**< The list of free surface edges (doublet of index in the nodesList. */

        std::vector<double> m_elementsDetJ;         /**< The Jacobian matrix determinant of each element. */
        std::vector<std::vector<std::vector<double>>> m_elementsInvJ;   /**< The inverse Jacobian matrix of each element. */
        std::vector<std::vector<std::vector<double>>> m_elementsJ;   /**< The Jacobian matrix of each element. */
        std::vector<double> m_freeSurfaceEdgeDetJ;         /**< The Jacobian matrix determinant of each free surface edge. */

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
        void computeElementsDetJ();

        /**
         * \brief Compute the inverse Jacobian matrix for gauss integration for each
         *        triangle.
         */
        void computeElementsInvJ();

        /**
         * \brief Compute the inverse Jacobian matrix for gauss integration for each
         *        triangle.
         */
        void computeElementsJ();

//        /**
//         * \brief Compute the determinant of the of the Jacobian matrix for gauss
//         *        integration for each edge.
//         */
//        void computeFreeSurfaceEdgeDetJ();

         /**
         * \brief Compute the mesh dimension from the .msh file.
         * \return The mesh dimension (2 or 3).
         */
        unsigned short computeMeshDim() const;

        /**
         * \brief Remesh the nodes in nodesList using CGAL (Delaunay triangulation
         *        and alpha-shape).
         */
        void triangulateAlphaShape();

        /**
         * \brief Remesh the nodes in nodesList using CGAL (Delaunay triangulation
         *        and alpha-shape) (2D).
         */
        void triangulateAlphaShape2D();

        /**
         * \brief Remesh the nodes in nodesList using CGAL (Delaunay triangulation
         *        and alpha-shape) (3D).
         */
        void triangulateAlphaShape3D();

        /**
         * \brief Removes nodes if they are too close from each other
         *       (d_nodes < gamma*hchar).
         * \return true if at least one node was deleted, false otherwise.
         */
        bool removeNodes();
};

#include "Mesh.inl"

#endif // MESH_HPP_INCLUDED
