#ifndef PARAMS_H_INCLUDED
#define PARAMS_H_INCLUDED

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

class Params
{
    public:
        Params();
        ~Params();

        void loadFromFile(std::string fileName);

        inline nlohmann::json getJSON() const { return m_json; }
;

    private:
        nlohmann::json m_json;
};

#endif // PARAMS_H_INCLUDED
