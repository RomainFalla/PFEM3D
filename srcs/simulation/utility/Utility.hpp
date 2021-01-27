#pragma once
#ifndef UTILITY_HPP_INCLUDED
#define UTILITY_HPP_INCLUDED

#include <iostream>
#include <bitset>

inline void printSize()
{
    std::cout << "Size of: " << "\n"
              << "\t - char: "              << sizeof(char) << "\n"
              << "\t - short: "             << sizeof(short) << "\n"
              << "\t - int: "               << sizeof(int) << "\n"
              << "\t - long: "              << sizeof(long) << "\n"
              << "\t - long long: "         << sizeof(long long) << "\n"
              << "\t - std::size_t: "       << sizeof(std::size_t) << "\n"
              << "\t - std::ptrdiff_t: "    << sizeof(std::ptrdiff_t) << "\n"
              << "\t - float: "             << sizeof(float) << "\n"
              << "\t - double: "            << sizeof(double) << "\n"
              << "\t - long double: "       << sizeof(long double) << "\n"
              << "\t - void*: "             << sizeof(void*) << "\n"
              << "\t - std::bitset<4>: "    << sizeof(std::bitset<4>) << "\n"
              << "\t - std::bitset<8>: "    << sizeof(std::bitset<8>) << std::endl;
}

#endif // UTILITY_HPP_INCLUDED

