#pragma once
#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include <string>
#include <map>

#include "Node.hpp"
#include "Element.hpp"
#include "Facet.hpp"

#include "mesh_defines.h"

struct MeshCreateInfo
{
    double hchar = 0;
    double alpha = 1;
    double gamma = 0;
    double omega = 1e16;
    std::vector<double> boundingBox = {};
    std::string mshFile = {};
};

/**
 * \class Mesh
 * \brief Represents a lagrangian mesh.
 */
class MESH_API Mesh
{
    public:
        Mesh()                              = delete;
        Mesh(const MeshCreateInfo& meshInfos);
        Mesh(const Mesh& mesh)              = delete;
        Mesh& operator=(const Mesh& mesh)   = delete;
        Mesh(Mesh&& mesh)                   = delete;
        Mesh& operator=(Mesh&& mesh)        = delete;
        ~Mesh()                             = default;

        /**
         * \brief Display the mesh parameters to console.
         */
        void displayToConsole() const noexcept;

        /**
         * \return The mesh dimension.
         */
        inline unsigned short getDim() const noexcept;

        /**
         * \param elm The index of the element.
         * \return A reference to the element.
         */
        inline const Element& getElement(std::size_t elm) const noexcept;

        /**
         * \return The number of elements in the mesh.
         */
        inline std::size_t getElementsCount() const noexcept;

        /**
         * \param face The index of the face.
         * \return A reference to the face.
         */
        inline const Facet& getFacet(std::size_t facet) const noexcept;

        /**
         * \return The number of boundary face in the mesh.
         */
        inline std::size_t getFacetsCount() const noexcept;

        /**
         * \param dimension The dimension of the reference element on which you want the weights.
         * \param n The requested number of Gauss points.
         * \return The Gauss points (x, y, z).
         */
        std::vector<std::array<double, 3>> getGaussPoints(uint8_t dimensiont, uint8_t n) const;

        /**
         * \param dimension The dimension of the reference element on which you want the weights.
         * \param n The requested number of Gauss points.
         * \return The weight associated to those Gauss points.
         */
        std::vector<double> getGaussWeight(uint8_t dimension, uint8_t n) const;

        /**
         * \return The characteristic size of the mesh.
         */
        inline double getHchar() const noexcept;

        /**
         * \return The number of nodes in the mesh.
         */
        inline const Node& getNode(std::size_t nodeIndex) const noexcept;

        /**
         * \return The number of nodes in the mesh.
         */
        inline std::size_t getNodesCount() const noexcept;

        /**
         * \param nodeIndex The index of the boundary node in the nodes list.
         * \return The initial position of the node.
         */
        inline std::array<double, 3> getBoundNodeInitPos(std::size_t nodeIndex) const;

        /**
         * \param nodeIndex The index of the boudary or free surface node in the nodes list.
         * \return The exterior normal of the boundary at the node.
         */
        inline std::array<double, 3> getFreeSurfaceNormal(std::size_t nodeIndex) const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \return The curvature of the boundary at the node.
         */
        inline double getFreeSurfaceCurvature(std::size_t nodeIndex) const;

        /**
         * \param nodeIndex The index of the node in the nodes list.
         * \return The physical group of that node.
         */
        inline std::string getNodeType(std::size_t nodeIndex) const noexcept;

        /**
         * \param dimension The dimension of the reference element on which you want the weights.
         * \return The size of the element in the reference coordinate system.
         */
        double getRefElementSize(uint8_t dimension) const;

        /**
         * \brief Perform remeshing on the mesh.
         */
        void remesh(bool verboseOutput);

        /**
         * \brief Restore the current nodes list from the list saved in saveNodesList.
         */
        void restoreNodesList();

        /**
         * \brief Save the current nodes list in another variable.
         */
        void saveNodesList();

        /**
         * \brief Activate or not the computation of normals and curvature (default is true).
         * \param activate should the computation be activated
         */
        inline void setComputeNormalCurvature(bool activate) noexcept;

        /**
         * \brief Set if the node is in a Dirichlet BC (i.e. a BC which we impose speed but do not move!).
         * \param nodeIndex The index of the node in the internal nodes list;
         * \param isDirichlet true if the node is in a Dirichlet BC, false otherwise;
         */
        inline void setNodeIsFixed(std::size_t nodeIndex, bool isFixed) noexcept;

        /**
         * \brief Set the position of a node.
         * \param nodeIndex The index of the node in the nodes list.
         * \param stateIndex The index of the state.
         * \param state The new value of the state.
         */
        inline void setNodeState(std::size_t nodeIndex, uint16_t stateIndex, double state) noexcept;

        /**
         * \brief Set the number of states to be stored at node level.
         * \param statesNumber The number of state per nodes.
         */
        inline void setStatesNumber(uint16_t statesNumber);

        /**
         * \brief Set the number of states to be stored at node level.
         * \param nodeIndex The index of the node which will be tagged.
         * \param tag The (positive) tag to set.
         */
        void setUserDefTag(std::size_t nodeIndex, int16_t tag);

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

    private:
        double m_hchar; /**< Characteristic size of an element (same as in .geo file). */
        double m_alpha; /**< Alpha parameter of the alpha-shape algorithm (triangles are discared if  r_circumcircle > alpha*hchar). */
        double m_omega; /**< Control the addition of node if a triangle is too big (a node is added if A_triangle > omege*hchar^2). */
        double m_gamma; /**< Control the deletetion of node if two are too close to each other (a node is deleted if d_nodes < gamma*hchar). */
        std::vector<double> m_boundingBox; /**< Box delimiting the zone of nodes existence (format: [xmin, ymin, xmax, ymax]). */

        bool m_computeNormalCurvature;

        unsigned short m_dim;   /**< The mesh dimension. */

        std::vector<Node> m_nodesList;      /**< List of nodes of the mesh. */
        std::vector<Node> m_nodesListSave;  /**< A copy of the nodes list (usefull for non-linear algorithm). */
        std::vector<Element> m_elementsList;    /**< The list of elements. */
        std::vector<Facet> m_facetsList;          /**< The list of boundary facets. */

        std::vector<std::string> m_tagNames; /**< The name of the tag of the nodes. */
        std::map<std::size_t, std::array<double, 3>> m_boundaryInitialPos;
        std::map<std::size_t, std::array<double, 3>> m_freeSurfaceNormal;
        std::map<std::size_t, double> m_freeSurfaceCurvature;

        /**
         * \brief Add nodes in element whose area is too big (A_tringle > omega*hchar^2.
         * \return true if at least one node was added, false otherwise).
         */
        bool addNodes(bool verboseOutput);

        /**
         * \brief Check if a node is outside the bounding box and deletes it if so.
         * \return true if at least one node was deleted, false otherwise.
         */
        bool checkBoundingBox(bool verboseOutput) noexcept;

         /**
         * \brief Compute the mesh dimension from the .msh file.
         * \return The mesh dimension (2 or 3).
         */
        void computeMeshDim();

         /**
         * \brief Compute the normal and curvature of each boundary and free surface.
         */
        void computeFSNormalCurvature();

        /**
         * \brief Compute the normal and curvature of each boundary and free surface (2D).
         */
        void computeFSNormalCurvature2D();

        /**
         * \brief Compute the normal and curvature of each boundary and free surface (3D).
         */
        void computeFSNormalCurvature3D();

        /**
         * \brief Load the nodes from a file using gmsh.
         * \param fileName The name of the .msh file.
         */
        void loadFromFile(const std::string& fileName);

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
        bool removeNodes(bool verboseOutput) noexcept;
};

#include "Mesh.inl"

#endif // MESH_HPP_INCLUDED
