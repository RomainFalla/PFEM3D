#include "Facet.hpp"

inline uint8_t Facet::getNodesCount() const noexcept
{
    return m_nodesIndexes.size();
}

inline std::size_t Facet::getNodeIndex(uint8_t node) const noexcept
{
    return m_nodesIndexes[node];
}

inline std::size_t Facet::getOutNodeIndex() const noexcept
{
    return m_outNodeIndex;
}

inline double Facet::getDetJ() const noexcept
{
    return m_detJ;
}

inline double Facet::getJ(uint8_t i, uint8_t j) const noexcept
{
    return m_J[i][j];
}

inline double Facet::getInvJ(uint8_t i, uint8_t j) const noexcept
{
    return m_invJ[i][j];
}

inline bool operator==(const Facet& a, const Facet& b) noexcept
{
    return std::equal(a.m_nodesIndexes.cbegin(), a.m_nodesIndexes.cend(), b.m_nodesIndexes.cbegin());
}
