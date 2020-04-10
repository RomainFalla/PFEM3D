#pragma once
#ifndef GMSHEXTRACTOR_HPP_INCLUDED
#define GMSHEXTRACTOR_HPP_INCLUDED

#include "Extractor.hpp"

class SOLVER_NO_EXPORT GMSHExtractor : public Extractor
{
    public:
        GMSHExtractor(const Solver& solver, const std::string& outFileName, double timeBetweenWriting,
                      const std::vector<std::string>& whatToWrite,
                      std::vector<std::string> whatCanBeWritten, std::string writeAs);

        ~GMSHExtractor();

        void update() override;

    private:
        std::vector<std::string> m_whatCanBeWritten;
        std::vector<bool> m_whatToWrite;
        std::string m_writeAs;
};

#endif // GMSHEXTRACTOR_HPP_INCLUDED
