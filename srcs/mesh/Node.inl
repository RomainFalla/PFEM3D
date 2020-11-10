#include "Node.hpp"

inline std::array<double, 3> Node::getPosition() const noexcept
{
    return m_position;
}

inline double Node::getPosition(uint8_t xyz) const noexcept
{
    return m_position[xyz];
}

inline double Node::getState(uint16_t state) const noexcept
{
    return m_states[state];
}

inline int16_t Node::getUserDefTag() const noexcept
{
    return m_userDefTag;
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
