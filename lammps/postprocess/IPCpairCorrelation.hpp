#include <vector>
#include <string>
#include <cmath>

#include "IPC.hpp"

class IPCisotropicPairCorrelationFunction {

public:
    IPCisotropicPairCorrelationFunction(int pSamplesBetweenZeroAndOne, Triad const& pSimulationBoxSide, int pNparticles);
    void accumulate(std::vector<IPC> particles);
    double print(std::string const& outputFileName);

private:
    IPCisotropicPairCorrelationFunction();
    int nParticles;
    Triad boxSide;

    double binsBetweenZeroAndOne;
    int paircorrelationTotalSamplings;

    std::vector<double> g;

    inline void relativePBC(double & x) { x -= std::lround(x); }

};
