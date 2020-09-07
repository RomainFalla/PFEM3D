#pragma once
#ifndef ELEMENT_HPP_INCLUDED
#define ELEMENT_HPP_INCLUDED

#include <vector>

/**
 * \typedef IndexType
 * \brief Type of the vector index for nodes and element (change it depending on
 *        maximum number of elements you want.
 */
typedef unsigned int IndexType;

/**
 * \struct Element
 * \brief Represents an element in PFEM method.
 */
struct Element
{
    std::vector<IndexType> nodesIndexes = {};        /**< Indexes of the nodes in the nodes list. */

    double detJ = 0;                                /**< Determinant of the Jacobian matrix. */
    std::vector<std::vector<double>> J = {};          /**< Jacobian matrix. */
    std::vector<std::vector<double>> invJ = {};       /**< Inverse Jacobian matrix. */
};

inline bool operator==(const Element& a, const Element& b) noexcept
{
    if(std::equal(a.nodesIndexes.cbegin(), a.nodesIndexes.cend(), b.nodesIndexes.cbegin()))
        return true;
    else
        return false;
}

#endif // ELEMENT_HPP_INCLUDED
