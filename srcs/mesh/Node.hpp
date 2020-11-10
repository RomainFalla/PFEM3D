#pragma once
#ifndef NODE_HPP_INCLUDED
#define NODE_HPP_INCLUDED

#include <array>
#include <cstdint>
#include <vector>

class Mesh;

/**
 * \class Node
 * \brief Represents a node in the PFEM.
 */
class Node
{
    public:
        Node()                              = default;
        Node(const Node& node)              = default;
        Node(Node&& node)                   = default;
        ~Node()                             = default;

        inline std::array<double, 3> getPosition() const noexcept;
        inline double getPosition(uint8_t xyz) const noexcept;
        inline double getState(uint16_t state) const noexcept;
        inline int16_t getUserDefTag() const noexcept;

        inline bool isFree() const noexcept;
        inline bool isBound() const noexcept;
        inline bool isOnFreeSurface() const noexcept;
        inline bool isFixed() const noexcept;

        friend inline bool operator==(const Node& a, const Node& b) noexcept;

        Node& operator=(const Node& node)   = default;
        Node& operator=(Node&& node)        = default;

    private:
        std::array<double, 3> m_position;           /**< Position of the node in 2D or 3D */

        std::vector<double> m_states;               /**< Variables defined at the node. Can be u,v,p or u,v,p,rho, etc. */

        std::vector<std::size_t> m_neighbourNodes;  /**< Indexes in the nodes list of the neighbour nodes */
        std::vector<std::size_t> m_elements;        /**< Index in the elements list of the elements which have this node */
        std::vector<std::size_t> m_facets;           /**< Index in the faces list of the faces which have this node */

        bool m_isBound = false;                     /**< Is the node a wall node */
        bool m_isOnFreeSurface = false;             /**< Is the node on the free surface */
        bool m_isFixed = false;                     /**< Is the node afixed (has a speed but does not move) ? */

        int16_t m_tag = -1;                         /**< Identify to which Physical group this node belongs to*/
        int16_t m_userDefTag = -1;                  /**< User Define Tag (can be used to tag nodes based on BC to applu*/

        friend class Mesh;
};

#include "Node.inl"

#endif // NODE_HPP_INCLUDED
