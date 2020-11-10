#include "Element.hpp"

inline uint8_t Element::getNodeCount() const noexcept
{
    return m_nodesIndexes.size();
}

inline std::size_t Element::getNodeIndex(uint8_t node) const noexcept
{
    return m_nodesIndexes[node];
}

inline double Element::getDetJ() const noexcept
{
    return m_detJ;
}

inline double Element::getJ(uint8_t i, uint8_t j) const noexcept
{
    return m_J[i][j];
}

inline double Element::getInvJ(uint8_t i, uint8_t j) const noexcept
{
    return m_invJ[i][j];
}

inline bool operator==(const Element& a, const Element& b) noexcept
{
    return std::equal(a.m_nodesIndexes.cbegin(), a.m_nodesIndexes.cend(), b.m_nodesIndexes.cbegin());
}
