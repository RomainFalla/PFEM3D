#pragma once
#ifndef POINTEXTRACTOR_HPP_INCLUDED
#define POINTEXTRACTOR_HPP_INCLUDED

#include <fstream>

#include "Extractor.hpp"

class SOLVER_NO_EXPORT PointExtractor : public Extractor
{
    public:
        PointExtractor(const Solver& solver, const std::string& outFileName, double timeBetweenWriting,
                       unsigned short stateToWrite, const std::vector<std::vector<double>>& points);
        ~PointExtractor();

        void update() override;

    private:
        std::ofstream m_outFile;
        unsigned short m_stateToWrite;
        std::vector<std::vector<double>> m_points;

        bool findElementIndex(const Mesh& mesh, IndexType& elementIndex, const std::vector<double>& point);
};

#endif // POINTEXTRACTOR_HPP_INCLUDED
