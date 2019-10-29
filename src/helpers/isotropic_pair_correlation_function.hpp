#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include "../IPC.hpp"

// All lengths has to be in natural units, only x[] works in units in which BoxSide=1.
class IsotropicPairCorrelationFunction {
// Computes the g(1,2) coefficients in the spherical harmonics expansion.
    // This computes g_ooo in the same way the isotropic g(r) is computed, and all the other
    // components like coeff*g_ooo*<spherical harmonics>. The theory can be found in the Hansen&McDonald.
    // The formula seems to come from Street&Tildesley, Proceedings of the Royal Society
    // of London. Series A, Mathematical and PhysicalSciences, Vol. 348, No. 1655 (Apr. 6, 1976), pp. 485-510
private:
    int nParticles;
    double simulationBoxSide;

    int binsBetweenZeroAndOne;
    double totalBinsInSimulationBoxSide;
    int paircorrelationTotalBins, paircorrelationTotalSamplings;

    std::vector<double> g;

    inline void relativePBC(double & x) { x -= std::lround(x); }

public:


    void initialize(int pSamplesBetweenZeroAndOne, double pSimulationBoxSide, int pNparticles) {
        simulationBoxSide = pSimulationBoxSide;
        nParticles = pNparticles;
        binsBetweenZeroAndOne = pSamplesBetweenZeroAndOne;
        totalBinsInSimulationBoxSide = binsBetweenZeroAndOne * simulationBoxSide;
        paircorrelationTotalBins = int( .5*std::sqrt(3.)*totalBinsInSimulationBoxSide ) + 1; // half a diagonal times the sampling in a unit
        paircorrelationTotalSamplings = 0;

        g.resize(paircorrelationTotalBins, 0.);
    }

    template <typename IPCtype>
    void compute(std::vector<IPCtype> particles) {
        // consistency check --- this should probably be done using enable_if or other template-specific helpers, but this works and it's much more readable
        static_assert(std::is_base_of<IPCbase, IPCtype>::value, "FATAL MISUSE OF cell_lists::compileLists(std::vector<IPCtype>) : IPCtype must be IPCbase or derived");
        ++paircorrelationTotalSamplings;
        for (int i=0; i<nParticles-1; ++i) {
            for (int k=i+1; k<nParticles; ++k) {
                double rik[3];
                for (int d: {0, 1, 2}) {
                    rik[d] = particles[i].ipcCenter.x[d] - particles[k].ipcCenter.x[d];
                    relativePBC(rik[d]);
                }
                const double r = std::sqrt(rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2]);
                const int R = (int) (r * totalBinsInSimulationBoxSide); // this works because rik is in [0:1) units
                g[R]  += 1.;
            }
        }
    }

    void compute(std::vector<IPE> particles) {
        ++paircorrelationTotalSamplings;
        for (int i=0; i<nParticles-1; ++i) {
            for (int k=i+1; k<nParticles; ++k) {
                double rik[3];
                for (int d: {0, 1, 2}) {
                    rik[d] = particles[i].cmPosition[d] - particles[k].cmPosition[d];
                    relativePBC(rik[d]);
                }
                const double r = std::sqrt(rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2]);
                const int R = (int) (r * totalBinsInSimulationBoxSide); // this works because rik is in [0:1) units
                g[R]  += 1.;
            }
        }
    }

    double print(std::string const& outputFileName) {
        std::ofstream outputFile(outputFileName);
        outputFile << std::scientific << std::setprecision(8);
        const double binSize = 1./binsBetweenZeroAndOne;
        const double normalization = (6.*std::pow(simulationBoxSide,3))/(paircorrelationTotalSamplings*4.*std::pow(nParticles,2)*M_PI*std::pow(binSize,3));
        double integral = 0.;
        for (int i=0; i<paircorrelationTotalBins; i++)
        {
            const double r = (i+.5)*binSize;
            const double shellVolume = pow(double(i+1),3)-pow(double(i),3);
            const double gScalingFactor = normalization/ shellVolume;
            g[i] = g[i]*gScalingFactor;

            outputFile << r << "\t" << g[i] << std::endl;

            integral += g[i]*shellVolume;
        }
        outputFile.close();

        const double density = nParticles/std::pow(simulationBoxSide,3);
        const double missingShellVolumeScaling = std::pow(binSize,3)*4.*M_PI/3.;
        return density*integral*missingShellVolumeScaling;
    }

};
