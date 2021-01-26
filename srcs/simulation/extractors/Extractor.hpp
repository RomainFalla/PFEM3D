#pragma once
#ifndef EXTRACTOR_HPP_INCLUDED
#define EXTRACTOR_HPP_INCLUDED

#include <string>

#include "../simulation_defines.h"

class Problem;

/**
 * \class Extractor
 * \brief Extractor interface.
 */
class SIMULATION_API Extractor
{
    public:
        Extractor()                                      = delete;
        /**
         * \param pProblem A pointer to the problem from which data will be extracted.
         * \param outFileName The file name in which data will be written.
         * \param timeBetweenWriting The simulation time between each write.
         */
        Extractor(Problem* pProblem, const std::string& outFileName, double timeBetweenWriting);
        Extractor(const Extractor& extractor)            = delete;
        Extractor& operator=(const Extractor& extractor) = delete;
        Extractor(Extractor&& extractor)                 = delete;
        Extractor& operator=(Extractor&& extractor)      = delete;
        virtual ~Extractor();

        /// Update the extractor state and write data if necessary
        virtual void update();

    protected:
        Problem* m_pProblem;
        std::string m_outFileName;

        double m_timeBetweenWriting;
        double m_nextWriteTrigger;
};

#endif // EXTRACTOR_HPP_INCLUDED
