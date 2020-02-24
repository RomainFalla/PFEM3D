#include <iostream>
#include <fstream>
#include <cerrno> //Change with something more c++ ;-)

#include "Params.hpp"

Params::Params()
{

}

Params::~Params()
{

}

void Params::loadFromFile(std::string fileName)
{
    std::ifstream paramFile(fileName);

    if(!paramFile.is_open())
    	throw std::runtime_error("The params file cannot be read!");

    paramFile >> m_json;
    paramFile.close();
}
