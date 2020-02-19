#ifndef PARAMS_H_INCLUDED
#define PARAMS_H_INCLUDED

#include <string>
#include <vector>

class Params
{
    public:
        double hchar;
        double alphaHchar;
        double omegaH2;
        double gammaH;

        std::vector<double> fluidParameters;

        double gravity;

        double maxTimeStep;
        double simuTime;
        double simuTimeToWrite;

        bool adaptDT;

        std::vector<double> sideParams;

        bool loadFromFile(std::string fileName);

        //std::vector<double> remeshingParameters;
};


#endif // PARAMS_H_INCLUDED
