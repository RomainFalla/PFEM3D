#include <array>

template<typename T>
constexpr std::array<std::array<T, 2>, 3> GP2Dpoints = {{{1/6, 1/6}, {2/3, 1/6}, {1/6, 2/3}}};

template<typename T>
constexpr std::array<T, 3> GP2Dweight = {1/3, 1/3, 1/3};
