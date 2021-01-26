#include "Mesh.hpp"

#include <array>
#include <iostream>

inline unsigned short Mesh::getDim() const noexcept
{
    return m_dim;
}

inline const Element& Mesh::getElement(std::size_t elm) const noexcept
{
    return m_elementsList[elm];
}

inline std::size_t Mesh::getElementsCount() const noexcept
{
    return m_elementsList.size();
}

inline const Facet& Mesh::getFacet(std::size_t facet) const noexcept
{
    return m_facetsList[facet];
}

inline std::size_t Mesh::getFacetsCount() const noexcept
{
    return m_facetsList.size();
}

inline std::string Mesh::getFacetType(std::size_t facetIndex) const noexcept
{
    std::size_t nodeIndex = m_facetsList[facetIndex].m_nodesIndexes[0];
    return m_tagNames[m_nodesList[nodeIndex].m_tag];
}

inline double Mesh::getHchar() const noexcept
{
    return m_hchar;
}

inline const Node& Mesh::getNode(std::size_t nodeIndex) const noexcept
{
    return m_nodesList[nodeIndex];
}

inline std::size_t Mesh::getNodesCount() const noexcept
{
    return m_nodesList.size();
}

std::array<double, 3> Mesh::getBoundNodeInitPos(std::size_t nodeIndex) const
{
    auto it = m_boundaryInitialPos.find(nodeIndex);

    if(it == m_boundaryInitialPos.end())
        throw std::runtime_error("initial position is only available for boundary nodes!");
    else
        return it->second;
}


inline std::array<double, 3> Mesh::getBoundFSNormal(std::size_t nodeIndex) const
{
    auto it = m_boundFSNormal.find(nodeIndex);

    if(it == m_boundFSNormal.end())
        throw std::runtime_error("normal is only available for free surface nodes!");
    else
        return it->second;
}

inline double Mesh::getFreeSurfaceCurvature(std::size_t nodeIndex) const
{
    auto it = m_freeSurfaceCurvature.find(nodeIndex);

    if(it == m_freeSurfaceCurvature.end())
        throw std::runtime_error("curvature is only available for free surface nodes!");
    else
        return it->second;
}

inline unsigned short Mesh::getNodesPerElm() const noexcept
{
    return m_dim + 1;
}

inline unsigned short Mesh::getNodesPerFacet() const noexcept
{
    return m_dim;
}

inline std::string Mesh::getNodeType(std::size_t nodeIndex) const noexcept
{
    return m_tagNames[m_nodesList[nodeIndex].m_tag];
}

inline bool Mesh::isNormalCurvComputed() const noexcept
{
    return m_computeNormalCurvature;
}

inline void Mesh::setComputeNormalCurvature(bool activate) noexcept
{
    m_computeNormalCurvature = activate;
}

inline void Mesh::setNodeFlag(std::size_t nodeIndex, unsigned short flag) noexcept
{
    m_nodesList[nodeIndex].m_userDefFlags.set(flag, 1);
}

inline void Mesh::setNodeIsFixed(std::size_t nodeIndex, bool isFixed) noexcept
{
    m_nodesList[nodeIndex].m_isFixed = isFixed;
}

inline void Mesh::setNodeState(std::size_t nodeIndex, unsigned int stateIndex, double state) noexcept
{
    m_nodesList[nodeIndex].m_states[stateIndex] = state;
}

void Mesh::setStatesNumber(unsigned int statesNumber)
{
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        m_nodesList[n].m_states.resize(statesNumber);
    }
}
