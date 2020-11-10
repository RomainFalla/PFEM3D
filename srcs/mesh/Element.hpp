#pragma once
#ifndef ELEMENT_HPP_INCLUDED
#define ELEMENT_HPP_INCLUDED

#include <array>
#include <cstdint>
#include <vector>

class Mesh;
class Node;

/**
 * \class Element
 * \brief Represents an element in the PFEM.
 */
class Element
{
    public:
        Element()                        = default;
        Element(const Element& element)  = default;
        Element(Element&& element)       = default;
        ~Element()                       = default;

        inline uint8_t getNodeCount() const noexcept;
        inline std::size_t getNodeIndex(uint8_t node) const noexcept;

        inline double getDetJ() const noexcept;
        inline double getJ(uint8_t i, uint8_t j) const noexcept;
        inline double getInvJ(uint8_t i, uint8_t j) const noexcept;

        friend inline bool operator==(const Element& a, const Element& b) noexcept;

        Element& operator=(const Element& element)   = default;
        Element& operator=(Element&& element)        = default;

    private:
        std::vector<std::size_t> m_nodesIndexes;      /**< Indexes of the nodes in the nodes list which compose this element. */

        double m_detJ;                                /**< Determinant of the Jacobian matrix of the element. */
        std::array<std::array<double, 3>, 3> m_J;     /**< Jacobian matrix of the element. */
        std::array<std::array<double, 3>, 3> m_invJ;  /**< Inverse Jacobian matrix of the element. */

        void computeJ(const std::vector<Node>& nodesList);
        void computeDetJ();
        void computeInvJ();

        friend class Mesh;
};

#include "Element.inl"

#endif // ELEMENT_HPP_INCLUDED
