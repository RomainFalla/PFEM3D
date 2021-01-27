#pragma once
#ifndef MINMAXEXTRACTOR_HPP_INCLUDED
#define MINMAXEXTRACTOR_HPP_INCLUDED

#include <fstream>

#include "Extractor.hpp"

class SIMULATION_API MinMaxExtractor : public Extractor
{
    public:
        MinMaxExtractor()                                                  = delete;
        MinMaxExtractor(Problem* pProblem, const std::string& outFileName, double timeBetweenWriting,
                        unsigned short coordinate, const std::string& minMax);
        MinMaxExtractor(const MinMaxExtractor& minmaxExtractor)            = delete;
        MinMaxExtractor& operator=(const MinMaxExtractor& minmaxExtractor) = delete;
        MinMaxExtractor(MinMaxExtractor&& minmaxExtractor)                 = delete;
        MinMaxExtractor& operator=(MinMaxExtractor&& minmaxExtractor)      = delete;
        ~MinMaxExtractor() override;

        void update() override;

    private:
        std::ofstream m_outFile;
        unsigned short m_coordinate;
        std::string m_minMax;
};

#endif // MINMAXEXTRACTOR_HPP_INCLUDED
