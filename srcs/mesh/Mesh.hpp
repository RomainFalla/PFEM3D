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
         * \return The mesh dimension.
         */
        inline unsigned short getDim() const;


        /**
         * \param elm The index of the element.
         * \return The vector containing the index of the nodes in nodesList
         *         making a certain element.
         */
        inline std::vector<IndexType> getElement(IndexType elm) const;

        /**
         * \param elm The element index.
         * \return The determinant of the of the Jacobian matrix
         */
        inline double getElementDetJ(IndexType elm) const;

        /**
         * \param elm The element index.
         * \param i The row index in the Jacobian matrix.
         * \param j The column index in the Jacobian matrix.
         * \return The ij element of the Jacobian matrix for a certain element.
         */
        inline double getElementInvJ(IndexType elm, unsigned short i, unsigned short j) const;

        /**
         * \return The number of elements in the mesh.
         */
        inline IndexType getElementsNumber() const;

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
         * \return The number of nodes in the mesh.
         */
        inline IndexType getNodesNumber() const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \return The vector initial position of the node.
         */
        inline std::vector<double> getNodeInitialPosition(IndexType nodeIndex) const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \return The vector position of the node.
         */
        inline std::vector<double> getNodePosition(IndexType nodeIndex) const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \param coordinate The wanted coordinate (0, 1, ...).
         * \return The coordinate of the node.
         */
        inline double getNodePosition(IndexType nodeIndex, unsigned short coordinate) const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \param state The wanted sate (0, 1, 2, 3, ...).
         * \return The value of the state for that node.
         */
        inline double getNodeState(IndexType nodeIndex, unsigned short state) const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \return The physical group of that node.
         */
        inline std::string getNodeType(IndexType nodeIndex) const;

        /**
         * \return The characteristic node adding parameter.
         */
        inline double getOmega() const;

        /**
         * \return The size of the element in the reference coordinate system.
         */
        inline double getRefElementSize() const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \return true if the node is free, false otherwise.
         */
        inline bool isNodeFree(IndexType nodeIndex) const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \return true if the node is on the boundary, false otherwise.
         */
        inline bool isNodeBound(IndexType nodeIndex) const;

         /**
         * \param nodeIndex The index of the node in the nodes list.
         * \return true if the node is on a Dirichlet boundary, false otherwise.
         */
        inline bool isNodeDirichlet(IndexType nodeIndex) const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \return true is the node is on the free surface, false otherwise.
         */
        inline bool isNodeOnFreeSurface(IndexType nodeIndex) const;

        /**
         * \brief Load the nodes from a file using gmsh.
         * \param fileName The name of the .msh file.
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
         * \brief Set if the node is in a Dirichlet BC (i.e. a BC which we impose speed but do not move!).
         * \param nodeIndex The index of the node in the internal nodes list;
         * \param isDirichlet true if the node is in a Dirichlet BC, false otherwise;
         */
        inline void setNodeIsDirichlet(IndexType nodeIndex, bool isDirichlet);

        /**
         * \brief Update the nodes position.
         * \param deltaPos The variation of the coordinate.
         */
        void updateNodesPosition(std::vector<double> deltaPos);

        /**
         * \brief Update the nodes position (from the saved nodeS List).
         * \param deltaPos The variation of the coordinate.
         */
        void updateNodesPositionFromSave(std::vector<double> deltaPos);

        /**
         * \brief Set the position of a node.
         * \param nodeIndex The index of the node in the nodes list.
         * \param stateIndex The index of the state.
         * \param state The new value of the state.
         */
        inline void setNodeState(IndexType nodeIndex, unsigned short stateIndex, double state);

        /**
         * \brief Set the number of states to be stored at node level.
         * \param statesNumber The number of state per nodes.
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
        std::vector<std::vector<IndexType>> m_elementsList;    /**< The list of element (triplet of index in the nodesList. */
        //std::vector<std::vector<IndexType>> m_freeSurfaceEdgesList;   /**< The list of free surface edges (doublet of index in the nodesList. */

        std::vector<std::string> m_tagNames; /**< The name of the tag of the nodes */

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
        void computeMeshDim();

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
