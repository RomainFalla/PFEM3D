#pragma once
#ifndef MASSEXTRACTOR_HPP_INCLUDED
#define MASSEXTRACTOR_HPP_INCLUDED

#include <fstream>

#include "Extractor.hpp"

class SOLVER_NO_EXPORT MassExtractor : public Extractor
{
    public:
        MassExtractor(const Solver& solver, const std::string& outFileName, double timeBetweenWriting);
        ~MassExtractor();

        void update() override;

    private:
        std::ofstream m_outFile;
};

#endif // MASSEXTRACTOR_HPP_INCLUDED
