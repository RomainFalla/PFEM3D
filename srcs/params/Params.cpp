#include <iostream>
#include <fstream>
#include <cerrno> //Change with something more c++ ;-)

#include "nlohmann/json.hpp"
#include "Params.hpp"

bool Params::loadFromFile(std::string fileName)
{
    std::ifstream paramFile(fileName);

    if(!paramFile.is_open())
    {
        std::cerr << "Something went wrong when trying to read the file "
                  << fileName << std::endl
                  << std::strerror(errno) << std::endl;

        paramFile.close();
        return false;
    }

    nlohmann::json j;
    paramFile >> j;
    paramFile.close();

    this->maxTimeStep = j["general"]["simulationMaxTimeStep"];
    this->simuTime = j["general"]["simulationTime"];
    this->simuTimeToWrite = j["general"]["simulationTimeToWrite"];

    this->hchar = j["remeshing"]["hchar"];
    this->alphaHchar = j["remeshing"]["alpha"];
    this->alphaHchar = this->alphaHchar*this->hchar;
    this->omegaH2 = j["remeshing"]["omega"];
    this->omegaH2 = this->omegaH2*this->hchar*this->hchar;
    this->gammaH = j["remeshing"]["gamma"];
    this->gammaH = this->gammaH*this->hchar;

    this->fluidParameters = j["physics"]["fluidParameters"].get<std::vector<double>>();
    this->gravity = j["physics"]["gravity"];

    return true;
}
