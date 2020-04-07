#pragma once
#ifndef EXTRACTOR_HPP_INCLUDED
#define EXTRACTOR_HPP_INCLUDED

#include <string>

#include "../../mesh/Mesh.hpp"

#include "Solver_export.h"

class SOLVER_NO_EXPORT Extractor
{
    public:
        Extractor(const std::string& outFileName, double timeBetweenWriting);
        ~Extractor();

        virtual void update(const Mesh& mesh, double currentTime, unsigned int currentStep);

    protected:
        std::string m_outFileName;

        double m_timeBetweenWriting;
        double m_nextWriteTrigger;
};

#endif // EXTRACTOR_HPP_INCLUDED
