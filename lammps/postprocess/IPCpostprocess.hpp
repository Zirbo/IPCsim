#ifndef __IPCPOSTPROCESS_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_HEADER_INCLUDED__

/*---------------------------------------------------------------------------------------
 * Postprocess for Inverse Patchy Colloids simulations with variable number of patches.
 *
 * asda
 *---------------------------------------------------------------------------------------*/
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "IPC.hpp"
#include "IPCpostprocessPotential.hpp"

class IPCpostprocess {
public:
    IPCpostprocess(std::string const& trajFilename, std::string const& inputFilename, const std::string &potDirName);
    void run();


private:
    IPCpostprocess();

    // input trajectory
    std::ifstream trajectoryFile;
    // state point
    int nIPCs;
    Triad boxSide;
    Ensemble ipcs;
    VectorOfTriads ipcOrientations;
    // geometry
    double ipcRadius, patchRadius, patchEccentricity, interactionRange;

    IPCpotential potential;

    void readFirstConfiguration();
    bool readNewConfiguration();
    void readIPCconfiguration();
    inline void relativePBC(double & x) {  x -= std::round(x);  }
    void computeOrientations();
    void computePerfectOrientations();
};

#endif //__IPCPOSTPROCESS_HEADER_INCLUDED__
