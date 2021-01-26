#pragma once
#ifndef POINTEXTRACTOR_HPP_INCLUDED
#define POINTEXTRACTOR_HPP_INCLUDED

#include <fstream>
#include <vector>

#include "Extractor.hpp"

class Mesh;

class SIMULATION_API PointExtractor : public Extractor
{
    public:
        PointExtractor()                                                = delete;
        PointExtractor(Problem* pProblem, const std::string& outFileName, double timeBetweenWriting,
                       std::string whatToWrite, const std::vector<std::vector<double>>& points);
        PointExtractor(const PointExtractor& pointExtractor)            = delete;
        PointExtractor& operator=(const PointExtractor& pointExtractor) = delete;
        PointExtractor(PointExtractor&& pointExtractor)                 = delete;
        PointExtractor& operator=(PointExtractor&& pointExtractor)      = delete;
        ~PointExtractor() override;

        void update() override;

    private:
        std::ofstream m_outFile;
        std::string m_whatToWrite;
        std::vector<std::vector<double>> m_points;

        bool findElementIndex(const Mesh& mesh, std::size_t& elementIndex, const std::vector<double>& point);
};

#endif // POINTEXTRACTOR_HPP_INCLUDED
