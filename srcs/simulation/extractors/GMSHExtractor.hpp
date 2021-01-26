#pragma once
#ifndef GMSHEXTRACTOR_HPP_INCLUDED
#define GMSHEXTRACTOR_HPP_INCLUDED

#include <vector>

#include "Extractor.hpp"

/**
 * \class GMSHExtractor
 * \brief Extractor to save data inside gmsh files.
 */
class SIMULATION_API GMSHExtractor : public Extractor
{
    public:
        GMSHExtractor()                                              = delete;
        /**
         * \param pProblem A pointer to the problem from which data will be extracted.
         * \param outFileName The file name in which data will be written.
         * \param timeBetweenWriting The simulation time between each write.
         * \param whatToWrite A vector containing the name of which data to write.
         * \param writeAs Hot to save the data. Can be "Nodes", "Elements" or "NodesElements".
         */
        GMSHExtractor(Problem* pProblem, const std::string& outFileName, double timeBetweenWriting,
                      const std::vector<std::string>& whatToWrite, std::string writeAs);
        GMSHExtractor(const GMSHExtractor& gmshExtractor)            = delete;
        GMSHExtractor& operator=(const GMSHExtractor& gmshExtractor) = delete;
        GMSHExtractor(GMSHExtractor&& gmshExtractor)                 = delete;
        GMSHExtractor& operator=(GMSHExtractor&& gmshExtractor)      = delete;
        ~GMSHExtractor() override;

        /// Update the extractor state and write data if necessary
        void update() override;

    private:
        std::vector<std::string> m_whatToWrite;
        std::vector<std::string> m_meshToWrite;
        std::string m_writeAs;

        static unsigned int m_initialized;
};

#endif // GMSHEXTRACTOR_HPP_INCLUDED
