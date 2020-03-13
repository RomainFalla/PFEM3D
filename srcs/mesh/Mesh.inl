#include "Mesh.hpp"

inline double Mesh::getDetJ(std::size_t elm) const
{
    assert(elm < m_elementList.size() && "elm should be between 0 and size - 1 !");
    assert(!m_detJ.empty());
    return m_detJ[elm];
}

inline std::vector<std::size_t> Mesh::getElement(std::size_t elm) const
{
    assert(elm < m_elementList.size() && "elm should be between 0 and size - 1 !");
    return m_elementList[elm];
}

inline std::size_t Mesh::getElementNumber() const
{
    return m_elementList.size();
}

inline double Mesh::getInvJ(std::size_t elm, unsigned short i, unsigned short j) const
{
    assert(elm < m_elementList.size() && "elm should be between 0 and size - 1 !");
    assert(!m_invJ.empty());
    assert(i < 2 && j < 2 && "i and j should be 0 or 1 in 2D");
    return m_invJ[elm][i][j];
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
