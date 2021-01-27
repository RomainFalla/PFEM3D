#pragma once
#ifndef STATESFROMQ_HPP_INCLUDED
#define STATESFROMQ_HPP_INCLUDED

#include <cassert>
#include "../../mesh/Mesh.hpp"
#include <Eigen/Dense>

inline void setNodesStatesfromQ(Mesh* pMesh, const Eigen::VectorXd& q, unsigned int beginState, unsigned int endState)
{
    assert(static_cast<std::size_t>(q.rows()) == (endState - beginState + 1)*pMesh->getNodesCount());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < pMesh->getNodesCount() ; ++n)
    {
        for (unsigned int s = beginState ; s <= endState ; ++s)
            pMesh->setNodeState(n, s, q((s - beginState)*pMesh->getNodesCount() + n));
    }
}

inline Eigen::VectorXd getQFromNodesStates(Mesh* pMesh, unsigned int beginState, unsigned int endState)
{
    Eigen::VectorXd q((endState - beginState + 1)*pMesh->getNodesCount());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < pMesh->getNodesCount() ; ++n)
    {
        const Node& node = pMesh->getNode(n);
        for (unsigned int s = beginState ; s <= endState ; ++s)
            q((s - beginState)*pMesh->getNodesCount() + n) = node.getState(s);
    }

    return q;
}

inline Eigen::VectorXd getElementState(Mesh* pMesh, const Element& element, unsigned int state)
{
    Eigen::VectorXd stateVec(pMesh->getNodesPerElm());

    if(pMesh->getDim() == 2)
    {
        stateVec << element.getNode(0).getState(state),
                    element.getNode(1).getState(state),
                    element.getNode(2).getState(state);
    }
    else
    {
        stateVec << element.getNode(0).getState(state),
                    element.getNode(1).getState(state),
                    element.getNode(2).getState(state),
                    element.getNode(3).getState(state);
    }

    return stateVec;
}

#endif
