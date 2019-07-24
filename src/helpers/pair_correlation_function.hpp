#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include "../IPC.hpp"

// All lengths has to be in natural units, only x[] works in units in which BoxSide=1.
class PairCorrelationFunction {
// Computes the g(1,2) coefficients in the spherical harmonics expansion.
    // This computes g_ooo in the same way the isotropic g(r) is computed, and all the other
    // components like coeff*g_ooo*<spherical harmonics>. The theory can be found in the Hansen&McDonald.
    // The formula seems to come from Street&Tildesley, Proceedings of the Royal Society
    // of London. Series A, Mathematical and PhysicalSciences, Vol. 348, No. 1655 (Apr. 6, 1976), pp. 485-510
private:
    bool computeAnisotropic, computeThirdOrder;
    int maxOrder;

    int nParticles;
    double simulationBoxSide;

    int binsBetweenZeroAndOne;
    double totalBinsInSimulationBoxSide;
    int paircorrelationTotalBins, paircorrelationTotalSamplings;

    std::vector<std::vector<double>> g;

    inline void absolutePBC(double & x) { x-=std::floor(x); }
    inline void relativePBC(double & x) { x-=std::lround(x); }

public:


    void initialize(int pSamplesBetweenZeroAndOne, double pSimulationBoxSide, int pNparticles,
                                             bool pComputeAnisotropic = false, bool pComputeThirdOrder = false) {
        computeThirdOrder = pComputeThirdOrder;
        computeAnisotropic = pComputeAnisotropic;
        if (computeThirdOrder)
            maxOrder = 30;
        else if (computeAnisotropic)
            maxOrder = 14;
        else
            maxOrder = 1;
        simulationBoxSide = pSimulationBoxSide;
        nParticles = pNparticles;
        binsBetweenZeroAndOne = pSamplesBetweenZeroAndOne;
        totalBinsInSimulationBoxSide = binsBetweenZeroAndOne * simulationBoxSide;
        paircorrelationTotalBins = int( .5*std::sqrt(3.)*totalBinsInSimulationBoxSide ) + 1; // half a diagonal times the sampling in a unit

        for(int i = 0; i < maxOrder ; ++i )
        {
            g[i].resize(paircorrelationTotalBins);
        }
    }

    template <typename IPCtype>
    void compute(std::vector<IPCtype> particles) {
        ++paircorrelationTotalSamplings;
        for (int i=0; i<nParticles-1; ++i) {
            for (int k=i+1; k<nParticles; ++k) {
                double rik[3];
                for (int d: {0, 1, 2}) {
                    rik[d] = particles[i].ipcCenter.x[d] - particles[k].ipcCenter.x[d];
                    relativePBC(rik[d]);
                }
                const double r = std::sqrt(rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2]);
                const int R = int( r * totalBinsInSimulationBoxSide); // this works because rik is in [0:1) units
                for (int d: {0, 1, 2})
                    rik[d] /= r;

                ///////////////////////calcolo delle armoniche sferiche
                g[0][R]  += 1.;                                                 //Y00*Y00
                /*
                if (computeAnisotropic) {
                    double cosQ1 = rik*w[i];                //cos theta i
                    double cosQ2 = rik*w[k];                //cos theta k
                    double cosQ1_2 = cosQ1*cosQ1;           //cos theta i quadro
                    double cosQ2_2 = cosQ2*cosQ2;           //cos theta k quadro
                    //il seno tra 0 e 3.14 è sempre positivo quindi non mi serve fare distinzioni
                    double senQ1_2 = 1.-cosQ1_2;
                    double senQ2_2 = 1.-cosQ2_2;
                    double senQ1 = sqrt(senQ1_2);
                    double senQ2 = sqrt(senQ2_2);
                    //come ottenere delta phi? Sottraggo alle orientazioni la componente parallela a r12;
                    //il prodotto scalare dei due vettorini rimasti (normalizzati) è cos delta phi
                    space::vec n1 = w[i] - rik*cosQ1;                n1 /= sqrt(n1*n1);
                    space::vec n2 = w[k] - rik*cosQ2;                n2 /= sqrt(n2*n2);
                    double cosp = n1*n2;
                    //polinomi di Legendre pronti all'uso
                    double P2Q1 = (3.*cosQ1_2-1.);
                    double P2Q2 = (3.*cosQ2_2-1.);
                    double cos2p = 2.*cosp*cosp-1.;
                    g[1][R]  += cosQ1;                                              //Y10*Y00
                    g[2][R]  += cosQ2;                                              //Y00*Y10
                    g[3][R]  += cosQ1*cosQ2;                                        //Y10*Y10
                    g[4][R]  += senQ1*senQ2*cosp;                                   //Y11*Y1-1
                    g[5][R]  += P2Q1;                                               //Y20*Y00
                    g[6][R]  += P2Q1*cosQ2;                                         //Y20*Y10
                    g[7][R]  += P2Q1*P2Q2;                                          //Y20*Y20
                    g[8][R]  += cosQ1*P2Q2;                                         //Y10*Y20
                    g[9][R]  += P2Q2;                                               //Y00*Y20
                    g[10][R] += cosQ1*senQ1*senQ2*cosp;                             //Y21*Y1-1
                    g[11][R] += cosQ1*cosQ2*senQ1*senQ2*cosp;                       //Y21*Y2-1
                    g[12][R] += cosQ2*senQ1*senQ2*cosp;                             //Y11*Y2-1
                    g[13][R] += senQ1_2*senQ2_2*cos2p;                              //Y22*Y2-2
                }
                if (computeThirdOrder) {
                    // compute Legendre polynomials
                    double P3Q1 = (5.*cosQ1_2-3.)*cosQ1;
                    double P3Q2 = (5.*cosQ2_2-3.)*cosQ2;
                    double cos3p = cosp*cos2p-sqrt((1.-cosp*cosp)*(1.-cos2p*cos2p));

                    g[14][R] += P3Q1;                                               //Y30*Y00
                    g[15][R] += P3Q1*cosQ2;                                         //Y30*Y10
                    g[16][R] += P3Q1*P2Q2;                                          //Y30*Y20
                    g[17][R] += P3Q1*P3Q2;                                          //Y30*Y30
                    g[18][R] += P2Q1*P3Q2;                                          //Y20*Y30
                    g[19][R] += cosQ1*P3Q2;                                         //Y10*Y30
                    g[20][R] += P3Q2;                                               //Y00*Y30

                    g[21][R] += senQ1*(5.*cosQ1_2-1.)*senQ2*cosp;                   //Y31*Y1-1
                    g[22][R] += senQ1*(5.*cosQ1_2-1.)*senQ2*cosQ2*cosp;             //Y31*Y2-1
                    g[23][R] += senQ1_2*cosQ1*senQ2_2*cos2p;                        //Y32*Y2-2
                    g[24][R] += senQ1*(5.*cosQ1_2-1.)*senQ2*(5.*cosQ2_2-1.)*cosp;   //Y31*Y3-1
                    g[25][R] += senQ1_2*cosQ1*senQ2_2*cosQ2*cos2p;                  //Y32*Y3-2
                    g[26][R] += senQ1_2*senQ1*senQ2_2*senQ2*cos3p;                  //Y33*Y3-3
                    g[27][R] += senQ1_2*senQ2_2*cosQ2*cos2p;                        //Y21*Y3-2
                    g[28][R] += senQ1*cosQ1*senQ2*(5.*cosQ2_2-1.)*cosp;             //Y22*Y3-1
                    g[29][R] += senQ1*senQ2*(5.*cosQ2_2-1.)*cosp;                   //Y11*Y3-1
                }*/
            ///////////////////////fine calcolo delle fottute armoniche sferiche
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
            g[0][i]        = g[0][i]*gScalingFactor; /*
            if (computeAnisotropic) {
                g[1][i]        = sqrt(3.)*g[1][i]*gScalingFactor*g[0][i];
                g[2][i]        = sqrt(3.)*g[2][i]*gScalingFactor*g[0][i];
                g[3][i]        = 3.*g[3][i]*gScalingFactor*g[0][i];
                g[4][i]        = -1.5*g[4][i]*gScalingFactor*g[0][i];
                g[5][i]        = 0.5*sqrt(5.)*g[5][i]*gScalingFactor*g[0][i];
                g[6][i]        = .5*sqrt(15)*g[6][i]*gScalingFactor*g[0][i];
                g[7][i]        = 1.25*g[7][i]*gScalingFactor*g[0][i];
                g[8][i]        = .5*sqrt(15)*g[8][i]*gScalingFactor*g[0][i];
                g[9][i]        = 0.5*sqrt(5.)*g[9][i]*gScalingFactor*g[0][i];
                g[10][i]       = -1.5*sqrt(5.)*g[10][i]*gScalingFactor*g[0][i];
                g[11][i]       = -7.5*g[11][i]*gScalingFactor*g[0][i];
                g[12][i]       = -1.5*sqrt(5.)*g[12][i]*gScalingFactor*g[0][i];
                g[13][i]       = 1.875*g[13][i]*gScalingFactor*g[0][i];
            }
            if(computeThirdOrder) {
                g[14][i]       = .5*sqrt(7.)*g[14][i]*gScalingFactor*g[0][i];
                g[15][i]       = .5*sqrt(21.)*g[15][i]*gScalingFactor*g[0][i];
                g[16][i]       = .25*sqrt(35.)*g[16][i]*gScalingFactor*g[0][i];
                g[17][i]       =  1.75*g[17][i]*gScalingFactor*g[0][i];
                g[18][i]       = .25*sqrt(35.)*g[18][i]*gScalingFactor*g[0][i];
                g[19][i]       = .5*sqrt(21.)*g[19][i]*gScalingFactor*g[0][i];
                g[20][i]       = .5*sqrt(7.)*g[20][i]*gScalingFactor*g[0][i];

                g[21][i]       = -.75*sqrt(3.5)*g[21][i]*gScalingFactor*g[0][i];
                g[22][i]       = -.75*sqrt(17.5)*g[22][i]*gScalingFactor*g[0][i];
                g[23][i]       = 1.875*sqrt(7.)*g[23][i]*gScalingFactor*g[0][i];
                g[24][i]       = -1.3125*g[24][i]*gScalingFactor*g[0][i];
                g[25][i]       = 13.125*g[25][i]*gScalingFactor*g[0][i];
                g[26][i]       = -2.1875*g[26][i]*gScalingFactor*g[0][i];
                g[27][i]       = 1.875*sqrt(7.)*g[27][i]*gScalingFactor*g[0][i];
                g[28][i]       = -.75*sqrt(17.5)*g[28][i]*gScalingFactor*g[0][i];
                g[29][i]       =  -.75*sqrt(3.5)*g[29][i]*gScalingFactor*g[0][i];
            }*/

            outputFile << r;
            for(int j = 0; j < maxOrder; ++j )
                outputFile << "\t" << g[j][i];
            outputFile << std::endl;

            integral += g[0][i]*shellVolume;
        }
        outputFile.close();

        return nParticles*std::pow(simulationBoxSide,3)*integral*4.*M_PI*std::pow(binSize,3)/3.;
    }

};
