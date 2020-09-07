#pragma once
#ifndef NODE_HPP_INCLUDED
#define NODE_HPP_INCLUDED

#include <vector>

#include "Element.hpp"

/**
 * \typedef IndexType
 * \brief Type of the vector index for nodes and element (change it depending on
 *        maximum number of elements you want.
 */
typedef unsigned int IndexType;

/**
 * \struct Node
 * \brief Represents a node in PFEM method.
 */
struct Node
{
    std::vector<double> position = {};          /**< Position of the node in 2D or 3D */

    std::vector<double> states = {};            /**< Variables defined at the node. Can be u,v,p or u,v,p,rho, etc. */

    std::vector<IndexType> neighbourNodes = {};       /**< Indexes in the nodes list of the neighbour nodes */
    std::vector<IndexType> belongingElements = {};    /**< Pointers to the elements owning that node */

    std::vector<double> initialPosition = {};   //Only for boundary nodes;

    bool isBound = false;                       /**< Is the node a wall node */
    bool isOnFreeSurface = false;               /**< Is the node on the free surface */
    bool isFree = false;                        /**< Is the node disconnected from any fluid elements */
    bool isDirichlet = false;                   /**< Is the node a Dirichlet BC node (has a speed but do not move) */

    unsigned short tag = 0;                         /**< Identify to which BC this node belongs to*/
};

inline bool operator==(const Node& a, const Node& b) noexcept
{
    if(std::equal(a.position.cbegin(), a.position.cend(), b.position.cbegin()))
        return true;
    else
        return false;
}

#endif // NODE_HPP_INCLUDED
