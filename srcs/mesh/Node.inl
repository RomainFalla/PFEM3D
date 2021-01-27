#include "Node.hpp"

inline std::array<double, 3> Node::getPosition() const noexcept
{
    return m_position;
}

inline double Node::getCoordinate(unsigned int xyz) const noexcept
{
    return m_position[xyz];
}

inline std::vector<double> Node::getStates() const noexcept
{
    return m_states;
}

inline double Node::getState(unsigned int state) const noexcept
{
    return m_states[state];
}

inline bool Node::getFlag(unsigned short flag) const noexcept
{
    return m_userDefFlags[flag];
}

inline unsigned int Node::getElementCount() const noexcept
{
    return m_elements.size();
}

inline unsigned int Node::getFacetCount() const noexcept
{
    return m_facets.size();
}

inline bool Node::isFree() const noexcept
{
    return m_elements.empty();
}

inline bool Node::isBound() const noexcept
{
    return m_isBound;
}

inline bool Node::isOnFreeSurface() const noexcept
{
    return m_isOnFreeSurface;
}

inline bool Node::isFixed() const noexcept
{
    return m_isFixed;
}

inline bool operator==(const Node& a, const Node& b) noexcept
{
    return std::equal(a.m_position.cbegin(), a.m_position.cend(), b.m_position.cbegin());
}
