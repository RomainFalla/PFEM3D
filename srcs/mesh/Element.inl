#include "Element.hpp"

inline double Element::getDetJ() const noexcept
{
    return m_detJ;
}

inline std::size_t Element::getNodeIndex(unsigned int node) const noexcept
{
    return m_nodesIndexes[node];
}

inline double Element::getJ(unsigned int i, unsigned int j) const noexcept
{
    return m_J[i][j];
}

inline double Element::getInvJ(unsigned int i, unsigned int j) const noexcept
{
    return m_invJ[i][j];
}

inline double Element::getLargestExtension() const noexcept 
{
    return m_largest_extension;
}

inline double Element::getCircumscribedRadius() const noexcept
{
    return m_r;
}

inline bool operator==(const Element& a, const Element& b) noexcept
{
    return std::equal(a.m_nodesIndexes.cbegin(), a.m_nodesIndexes.cend(), b.m_nodesIndexes.cbegin());
}


