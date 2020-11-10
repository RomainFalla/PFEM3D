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


inline std::array<double, 3> Mesh::getFreeSurfaceNormal(std::size_t nodeIndex) const
{
    auto it = m_freeSurfaceNormal.find(nodeIndex);

    if(it == m_freeSurfaceNormal.end())
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

inline std::string Mesh::getNodeType(std::size_t nodeIndex) const noexcept
{
    return m_tagNames[m_nodesList[nodeIndex].m_tag];
}

inline void Mesh::setComputeNormalCurvature(bool activate) noexcept
{
    m_computeNormalCurvature = activate;
}

inline void Mesh::setNodeIsFixed(std::size_t nodeIndex, bool isFixed) noexcept
{
    m_nodesList[nodeIndex].m_isFixed = isFixed;
}

inline void Mesh::setNodeState(std::size_t nodeIndex, uint16_t stateIndex, double state) noexcept
{
    m_nodesList[nodeIndex].m_states[stateIndex] = state;
}

void Mesh::setStatesNumber(unsigned short statesNumber)
{
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        m_nodesList[n].m_states.resize(statesNumber);
    }
}
