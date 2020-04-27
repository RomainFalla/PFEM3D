#ifndef NODE_HPP_INCLUDED
#define NODE_HPP_INCLUDED

#include <vector>

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
    std::vector<double> position;       /**< Position of the node in 2D or 3D */

    std::vector<double> states;         /**< Variables defined at the node. Can be u,v,p or u,v,p,rho, etc. */

    std::vector<IndexType> neighbourNodes;    /**< Indexes in the nodes list of the neighbour nodes */

    std::vector<double> initialPosition; //Only for boundary nodes;

    bool isBound;                       /**< Is the node a wall node */
    bool isOnFreeSurface;               /**< Is the node on the free surface */
    bool isFree;                        /**< Is the node disconnected from any fluid
                                             elements */
    bool isDirichlet;                   /**< Is the node a Dirichlet BC node (has a speed but do not move) */

    unsigned short tag;                   /**< Identify to which BC this node belongs to*/
};

inline bool operator==(const Node& a, const Node& b)
{
    for(unsigned short i = 0 ; i < a.position.size() ; ++i)
    {
        if(a.position[i] != b.position[i])
            return false;
    }

    return true;
}

#endif // NODE_HPP_INCLUDED
