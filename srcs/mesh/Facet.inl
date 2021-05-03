#include "Facet.hpp"

inline std::size_t Facet::getNodeIndex(unsigned int node) const noexcept
{
    return m_nodesIndexes[node];
}

inline std::size_t Facet::getOutNodeIndex() const noexcept
{
    return m_outNodeIndexes[0];
}

inline double Facet::getDetJ() const noexcept
{
    return m_detJ;
}

inline double Facet::getJ(unsigned int i, unsigned int j) const noexcept
{
    return m_J[i][j];
}

inline double Facet::getInvJ(unsigned int i, unsigned int j) const noexcept
{
    return m_invJ[i][j];
}

inline bool operator==(const Facet& a, const Facet& b) noexcept
{
    return std::equal(a.m_nodesIndexes.cbegin(), a.m_nodesIndexes.cend(), b.m_nodesIndexes.cbegin());
}
