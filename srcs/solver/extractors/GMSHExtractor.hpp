#pragma once
#ifndef GMSHEXTRACTOR_HPP_INCLUDED
#define GMSHEXTRACTOR_HPP_INCLUDED

#include "Extractor.hpp"

class SOLVER_NO_EXPORT GMSHExtractor : public Extractor
{
    public:
        GMSHExtractor(const std::string& outFileName, double timeBetweenWriting,
                      const std::vector<std::string>& whatToWrite,
                      std::vector<std::string> whatCanBeWritten, std::string writeAs,
                      unsigned short statesNumber);

        ~GMSHExtractor();

        void update(const Mesh& mesh, double currentTime, unsigned int currentStep) override;

    private:
        std::vector<std::string> m_whatCanBeWritten;
        std::vector<bool> m_whatToWrite;
        std::string m_writeAs;
        unsigned short m_statesNumber;
};

#endif // GMSHEXTRACTOR_HPP_INCLUDED
