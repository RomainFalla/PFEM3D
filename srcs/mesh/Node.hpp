#ifndef NODE_HPP_INCLUDED
#define NODE_HPP_INCLUDED

#include <cassert>
#include <vector>


/**
 * \struct Node
 * \brief Represents a node in PFEM method.
 */
struct Node
{
    std::vector<double> position;       /**< Position of the node in 2D or 3D*/

    std::vector<double> states;         /**< Variables defined at the node.
                                             Can be u,v,p or u,v,p,rho, etc.*/

    std::vector<std::size_t> neighbourNodes;    /**< Indexes in the nodes list of the
                                                     neighbour nodes*/

    bool isBound;                       /**< Is the node a wall node*/
    bool isOnFreeSurface;               /**< Is the node on the free surface*/
    bool isFree;                        /**< Is the node disconnected from any fluid
                                             elements*/

    bool touched;                       /**< Is the node next to a node which should
                                             be deleted */
    bool toBeDeleted;                   /**< Should the node be deleted */

    /**
     * \brief Constructor
     * \param dim The spatial dimension (2 or 3).
     * \param nUnknowns The number of unknowns of the problem.
     */
    Node(unsigned int dim, unsigned int nUnknowns)
    {
        assert(dim == 2 || dim == 3);
        this->position.resize(dim);
        this->states.resize(nUnknowns);

        this->isBound = false;
        this->isOnFreeSurface = false;
        this->isFree = true;

        this->touched = false;
        this->toBeDeleted = false;
    }
};

#endif // NODE_HPP_INCLUDED
