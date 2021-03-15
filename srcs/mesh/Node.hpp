#pragma once
#ifndef NODE_HPP_INCLUDED
#define NODE_HPP_INCLUDED

#include <array>
#include <bitset>
#include <cstdint>
#include <vector>

#include "mesh_defines.h"

class Mesh;
class Facet;
class Element;

/**
 * \class Node
 * \brief Represents a node in the PFEM.
 */
class MESH_API Node
{
    public:
        /// \param mesh a reference to the mesh
        Node(Mesh& mesh);
        Node(const Node& node)              = default;
        Node(Node&& node)                   = default;
        ~Node()                             = default;

        /// \return The distance between the two nodes n0 and n1
        static double distance(const Node& n0, const Node& n1);

        /// \return The (x, y, z) position of the node.
        inline std::array<double, 3> getPosition() const noexcept;

        /// \param xyz The index of the coordinate (x, y, z).
        /// \return The asked coordinate of the node.
        inline double getCoordinate(unsigned int xyz) const noexcept;

        /// \param state The index of the state.
        /// \return The asked state of the node.
        inline double getState(unsigned int state) const noexcept;

        /// \return The states of the node.
        inline std::vector<double> getStates() const noexcept;

        /// \return The Local mesh size at the node.
        inline double getLocalMeshSize() const noexcept;

        /// \set The Local mesh size = tmsize, at the node.
        inline void setLocalMeshSize(double tmsize);

        /// \return The user-defined flag of the node.
        inline bool getFlag(unsigned short flag) const noexcept;

        /// \param elementIndex The index of the node's element.
        /// \return A constant reference to the element.
        const Element& getElement(unsigned int elementIndex) const noexcept;

        /// \return The number of elements the node is inside.
        inline unsigned int getElementCount() const noexcept;

        /// \param facetIndex The index of the node's facet.
        /// \return A constant reference to the facet.
        const Facet& getFacet(unsigned int facetIndex) const noexcept;

        /// \return The number of facets the node is inside.
        inline unsigned int getFacetCount() const noexcept;

        /// \return Is the node on a boundary ?
        inline bool isBound() const noexcept;

        /// \return Is the node near a boundary ?
        bool isContact() const noexcept;

        /// \return Is the node position fixed ?
        inline bool isFixed() const noexcept;

        /// \return Is the node on any type of boundary ?
        inline bool isOnBoundary() const noexcept;

        /// \return Is the node free of elements ?
        inline bool isFree() const noexcept;

        /// \return Is the node on the free surface ?
        inline bool isOnFreeSurface() const noexcept;

        /// \return the natural mesh size at the node
        inline double getNaturalMeshSize() const noexcept;

        /// \update the natural mesh size at the node
        void updateNaturalMeshSize();



        /// \return the point weight used for the weighted triangulation and alpha shape
        double getWeight(double alphaRatio,double minTargetMeshSize);

        friend inline bool operator==(const Node& a, const Node& b) noexcept;

        Node& operator=(const Node& node)   = default;
        Node& operator=(Node&& node)        = default;

    private:
        Mesh* m_pMesh;                              /**< A pointer to the mesh from which the facet comes from. */

        std::array<double, 3> m_position;           /**< Position of the node in 2D or 3D. */

        std::vector<double> m_states;               /**< Variables defined at the node. Can be u,v,p or u,v,p,rho, etc. */

        std::vector<std::size_t> m_neighbourNodes;  /**< Indexes in the nodes list of the neighbour nodes. */
        std::vector<std::size_t> m_elements;        /**< Index in the elements list of the elements which have this node.*/
        std::vector<std::size_t> m_facets;           /**< Index in the faces list of the faces which have this node. */

        bool m_isBound = false;                     /**< Is the node a wall node. */
        bool m_isOnFreeSurface = false;             /**< Is the node on the free surface. */
        bool m_isFixed = false;                     /**< Is the node fixed (has a speed but does not move) ? */
        bool m_isOnBoundary = false;                /**< Is the node of any boundary type. */

        double m_target_mesh_size = 0.;
        double m_natural_mesh_size = 0.;

        int m_tag = -1;                             /**< Identify to which Physical group this node belongs to.*/
        std::bitset<8> m_userDefFlags;              /**< User define flags (can be used to tag nodes based on BC to apply).*/

        friend class Mesh;
        friend class Element;
};

#include "Node.inl"

#endif // NODE_HPP_INCLUDED
