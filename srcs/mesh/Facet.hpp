#pragma once
#ifndef FACET_HPP_INCLUDED
#define FACET_HPP_INCLUDED

#include <array>
#include <cstdint>
#include <vector>

class Mesh;
class Node;

/**
 * \class Facet
 * \brief Represents an Facet in the PFEM.
 */
class Facet
{
    public:
        Facet()                        = default;
        Facet(const Facet& facet)      = default;
        Facet(Facet&& facet)           = default;
        ~Facet()                       = default;

        inline uint8_t getNodesCount() const noexcept;
        inline std::size_t getNodeIndex(uint8_t node) const noexcept;
        inline std::size_t getOutNodeIndex() const noexcept;

        inline double getDetJ() const noexcept;
        inline double getJ(uint8_t i, uint8_t j) const noexcept;
        inline double getInvJ(uint8_t i, uint8_t j) const noexcept;

        friend inline bool operator==(const Facet& a, const Facet& b) noexcept;

        Facet& operator=(const Facet& facet)   = default;
        Facet& operator=(Facet&& facet)        = default;

    private:
        std::vector<std::size_t> m_nodesIndexes;      /**< Indexes of the nodes in the nodes list which compose this facet. */
        std::size_t m_outNodeIndex;                   /**< Indexes of the node which is "in front of" this facet. */

        double m_detJ;                                /**< Determinant of the Jacobian matrix of the facet. */
        std::array<std::array<double, 2>, 3> m_J;     /**< Jacobian matrix of the facet. */
        std::array<std::array<double, 3>, 2> m_invJ;  /**< Inverse Jacobian matrix of the facet. */

        void computeJ(const std::vector<Node>& nodesList);
        void computeDetJ();
        void computeInvJ(const std::vector<Node>& nodesList);

        friend class Mesh;
};

#include "Facet.inl"

#endif // FACET_HPP_INCLUDED
