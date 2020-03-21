#include "Mesh.hpp"

inline double Mesh::getAlpha() const
{
    return m_p.alpha;
}

inline std::vector<std::size_t> Mesh::getElement(std::size_t elm) const
{
    assert(elm < m_elementsList.size() && "elm should be between 0 and size - 1 !");
    return m_elementsList[elm];
}

inline double Mesh::getElementDetJ(std::size_t elm) const
{
    assert(elm < m_elementsList.size() && "elm should be between 0 and size - 1 !");
    assert(!m_elementDetJ.empty());
    return m_elementDetJ[elm];
}

inline double Mesh::getElementInvJ(std::size_t elm, unsigned short i, unsigned short j) const
{
    assert(elm < m_elementsList.size() && "elm should be between 0 and size - 1 !");
    assert(!m_elementInvJ.empty());
    assert(i < 2 && j < 2 && "i and j should be 0 or 1 in 2D");
    return m_elementInvJ[elm][i][j];
}

inline std::size_t Mesh::getElementsNumber() const
{
    return m_elementsList.size();
}

inline double Mesh::getFreeSurfaceDetJ(std::size_t edge) const
{
    assert(edge < m_freeSurfaceEdgesList.size() && "edge should be between 0 and size - 1 !");
    return m_freeSurfaceEdgeDetJ[edge];


}

inline std::vector<std::size_t>  Mesh::getFreeSurfaceEdge(std::size_t edge) const
{
    return m_freeSurfaceEdgesList[edge];
}

inline std::size_t  Mesh::getFreeSurfaceEdgesNumber() const
{
    return m_freeSurfaceEdgesList.size();
}

inline double Mesh::getGamma() const
{
    return m_p.gamma;
}

inline double Mesh::getHchar() const
{
    return m_p.hchar;
}

inline unsigned short Mesh::getMeshDim() const
{
    return m_dim;
}

inline std::size_t Mesh::getNodesNumber() const
{
    return m_nodesList.size();
}

inline double Mesh::getNodePosition(std::size_t nodeIndex, unsigned short coordinate) const
{
    return m_nodesList[nodeIndex].position[coordinate];
}

inline double Mesh::getNodeState(std::size_t nodeIndex, unsigned short state) const
{
    return m_nodesList[nodeIndex].states[state];
}

inline double Mesh::getOmega() const
{
    return m_p.omega;
}

inline bool Mesh::isNodeFree(std::size_t nodeIndex) const
{
    return m_nodesList[nodeIndex].isFree;
}

inline bool Mesh::isNodeBound(std::size_t nodeIndex) const
{
    return m_nodesList[nodeIndex].isBound;
}

inline bool Mesh::isNodeFluidInput(std::size_t nodeIndex) const
{
    return m_nodesList[nodeIndex].isFluidInput;
}

inline void Mesh::setNodeState(std::size_t nodeIndex, unsigned short stateIndex, double state)
{
    m_nodesList[nodeIndex].states[stateIndex] = state;
}
