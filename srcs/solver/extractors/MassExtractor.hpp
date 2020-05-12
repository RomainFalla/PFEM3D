#pragma once
#ifndef MASSEXTRACTOR_HPP_INCLUDED
#define MASSEXTRACTOR_HPP_INCLUDED

#include <fstream>

#include "Extractor.hpp"

class SOLVER_API MassExtractor : public Extractor
{
    public:
        MassExtractor()                                              = delete;
        MassExtractor(const Solver& solver, const std::string& outFileName, double timeBetweenWriting);
        MassExtractor(const MassExtractor& massExtractor)            = delete;
        MassExtractor& operator=(const MassExtractor& massExtractor) = delete;
        MassExtractor(MassExtractor&& massExtractor)                 = delete;
        MassExtractor& operator=(MassExtractor&& massExtractor)      = delete;
        ~MassExtractor() override;

        void update() override;

    private:
        std::ofstream m_outFile;
};

#endif // MASSEXTRACTOR_HPP_INCLUDED
