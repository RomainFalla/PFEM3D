#pragma once
#ifndef CLOCK_HPP_INCLUDED
#define CLOCK_HPP_INCLUDED

#include <chrono>
#include <iostream>

class Clock
{
    public:
        Clock() = default;
        ~Clock() = default;

        void start()
        {
            m_start = std::chrono::high_resolution_clock::now();
        }

        void end()
        {
            m_end = std::chrono::high_resolution_clock::now();
        }

        void displayDT(const std::string& text)
        {
            auto ellapsedTimeMeasure = std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start);
            std::cout << text << static_cast<double>(ellapsedTimeMeasure.count())/1000.0 << " s" << std::endl;
        }

    private:
        std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
        std::chrono::time_point<std::chrono::high_resolution_clock> m_end;
};

#endif // CLOCK_HPP_INCLUDED

