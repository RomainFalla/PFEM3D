#pragma once
#ifndef MINMAXEXTRACTOR_HPP_INCLUDED
#define MINMAXEXTRACTOR_HPP_INCLUDED

#include <fstream>

#include "Extractor.hpp"

class SOLVER_NO_EXPORT MinMaxExtractor : public Extractor
{
    public:
        MinMaxExtractor(const std::string& outFileName, double timeBetweenWriting,
                        unsigned short coordinate, const std::string& minMax, unsigned short meshDim);
        ~MinMaxExtractor();

        void update(const Mesh& mesh, double currentTime, unsigned int currentStep) override;

    private:
        std::ofstream m_outFile;
        unsigned short m_coordinate;
        std::string m_minMax;
};

#endif // MINMAXEXTRACTOR_HPP_INCLUDED
