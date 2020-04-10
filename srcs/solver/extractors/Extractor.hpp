#pragma once
#ifndef EXTRACTOR_HPP_INCLUDED
#define EXTRACTOR_HPP_INCLUDED

#include <string>

#include "../../mesh/Mesh.hpp"

#include "Solver_export.h"

class Solver;

class SOLVER_NO_EXPORT Extractor
{
    public:
        Extractor(const Solver& solver, const std::string& outFileName, double timeBetweenWriting);
        ~Extractor();

        virtual void update();

    protected:
        const Solver& m_solver;
        std::string m_outFileName;

        double m_timeBetweenWriting;
        double m_nextWriteTrigger;
};

#endif // EXTRACTOR_HPP_INCLUDED
