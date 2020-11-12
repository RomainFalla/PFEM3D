#pragma once
#ifndef EXTRACTOR_HPP_INCLUDED
#define EXTRACTOR_HPP_INCLUDED

#include <string>

#include "../solver_defines.h"

class Solver;

class SOLVER_API Extractor
{
    public:
        Extractor()                                      = delete;
        Extractor(const Solver& solver, const std::string& outFileName, double timeBetweenWriting);
        Extractor(const Extractor& extractor)            = delete;
        Extractor& operator=(const Extractor& extractor) = delete;
        Extractor(Extractor&& extractor)                 = delete;
        Extractor& operator=(Extractor&& extractor)      = delete;
        virtual ~Extractor();

        virtual void update();

    protected:
        const Solver& m_solver;
        std::string m_outFileName;

        double m_timeBetweenWriting;
        double m_nextWriteTrigger;
};

#endif // EXTRACTOR_HPP_INCLUDED
