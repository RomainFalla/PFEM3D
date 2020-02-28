#ifndef PARAMS_H_INCLUDED
#define PARAMS_H_INCLUDED

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

/**
 * \class Params
 * \brief Stores the json object containing the params (as well as the number of
 *        OMP threads).
 */
class Params
{
    public:
        Params();
        ~Params();

        /**
         * \brief Load in an object the parameters from a json file.
         * \param fileName The name of the parameters json file.
         */
        void loadFromFile(std::string fileName);

        /**
         * \brief Return an object containing the parameters.
         * \return The nlohmann::json object contaning the parameters.
         */
        inline nlohmann::json getJSON() const { return m_json; }

        /**
         * \brief Return the number of OpenMP threads the program is using.
         * \return The number of OpenMP threads the program is using.
         */
        inline unsigned int getNumOMPThreads() const { return m_threadsOMP; }

    private:
        nlohmann::json m_json; /**< A json object containing the parameters. */
        unsigned int m_threadsOMP; /**< The number of OpenMP threads the program use. */
};

#endif // PARAMS_H_INCLUDED
