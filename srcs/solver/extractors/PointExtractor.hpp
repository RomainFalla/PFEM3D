#pragma once
#ifndef POINTEXTRACTOR_HPP_INCLUDED
#define POINTEXTRACTOR_HPP_INCLUDED

#include <fstream>
#include <vector>

#include "Extractor.hpp"

class Mesh;

class SOLVER_API PointExtractor : public Extractor
{
    public:
        PointExtractor()                                                = delete;
        PointExtractor(const Solver& solver, const std::string& outFileName, double timeBetweenWriting,
                       unsigned short stateToWrite, const std::vector<std::vector<double>>& points);
        PointExtractor(const PointExtractor& pointExtractor)            = delete;
        PointExtractor& operator=(const PointExtractor& pointExtractor) = delete;
        PointExtractor(PointExtractor&& pointExtractor)                 = delete;
        PointExtractor& operator=(PointExtractor&& pointExtractor)      = delete;
        ~PointExtractor() override;

        void update() override;

    private:
        std::ofstream m_outFile;
        unsigned short m_stateToWrite;
        std::vector<std::vector<double>> m_points;

        bool findElementIndex(const Mesh& mesh, std::size_t& elementIndex, const std::vector<double>& point);
};

#endif // POINTEXTRACTOR_HPP_INCLUDED
