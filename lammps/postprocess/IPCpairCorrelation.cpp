#include <fstream>
#include <iomanip>
#include <cmath>
#include "IPCpairCorrelation.hpp"

IPCisotropicPairCorrelationFunction::IPCisotropicPairCorrelationFunction(int pSamplesBetweenZeroAndOne, Triad const& pSimulationBoxSide, int pNparticles) {
    boxSide = pSimulationBoxSide;
    nParticles = pNparticles;
    binsBetweenZeroAndOne = pSamplesBetweenZeroAndOne;
    paircorrelationTotalSamplings = 0;

    const double boxDiagonal = std::sqrt( std::pow(boxSide[0],2) + std::pow(boxSide[1],2) + std::pow(boxSide[2],2) );
    const int pairCorrelationTotalBins = int( binsBetweenZeroAndOne*0.5*boxDiagonal ) + 1;
    g.resize(pairCorrelationTotalBins, 0.);
}

void IPCisotropicPairCorrelationFunction::accumulate(std::vector<IPC> particles) {
    ++paircorrelationTotalSamplings;
    for (int i=0; i<nParticles-1; ++i) {
        for (int j=i+1; j<nParticles; ++j) {
            double r = 0.0;
            for (int d: {0, 1, 2}) {
                double rij = particles[i].ipcCenter.x[d] - particles[j].ipcCenter.x[d];
                relativePBC(rij);
                r += std::pow(rij*boxSide[d], 2);
            }
            r = std::sqrt(r);
            const int R = (int) (r * binsBetweenZeroAndOne);
            g[R]  += 1.;
        }
    }
}

double IPCisotropicPairCorrelationFunction::print(std::string const& outputFileName) {
    std::ofstream outputFile(outputFileName);
    outputFile << std::scientific << std::setprecision(8);
    const double binSize = 1./binsBetweenZeroAndOne;
    const double boxVolume = boxSide[0]*boxSide[1]*boxSide[2];
    const double normalization = (6.*boxVolume)/(paircorrelationTotalSamplings*4.*std::pow(nParticles,2)*M_PI*std::pow(binSize,3));
    double integral = 0.;
    for (int i = 0; i < int(g.size()); ++i)
    {
        const double r = (i+.5)*binSize;
        const double shellVolume = pow(double(i+1),3)-pow(double(i),3);
        const double gScalingFactor = normalization/ shellVolume;
        g[i] = g[i]*gScalingFactor;

        outputFile << r << "\t" << g[i] << std::endl;

        integral += g[i]*shellVolume;
    }
    outputFile.close();

    const double density = nParticles/boxVolume;
    const double missingShellVolumeScaling = std::pow(binSize,3)*4.*M_PI/3.;
    return density*integral*missingShellVolumeScaling;
}
