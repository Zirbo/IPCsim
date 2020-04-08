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
    double boxSideX, boxSideY, boxSideZ;
    std::vector<IPC> particles;
    // geometry
    double ipcRadius, patchRadius, patchEccentricity, interactionRange;

    IPCpotential potential;

    void readFirstConfiguration();
    bool readNewConfiguration();
    void readIPCconfiguration();
};

#endif //__IPCPOSTPROCESS_HEADER_INCLUDED__
