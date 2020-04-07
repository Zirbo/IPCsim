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


class IPCpostprocess {
public:
    IPCpostprocess(std::string const& trajFilename, std::string const& inputFilename);
    void run();


private:
    IPCpostprocess();

    // input trajectory
    std::ifstream trajectoryFile;
    // state point
    int nIPCs;
    double boxSideX, boxSideY, boxSideZ;
    // geometry
    double ipcRadius, patchRadius, patchEccentricity, interactionRange;


    std::vector<IPC> particles;

    void readFirstConfiguration();
    bool readNewConfiguration();
    void readIPCconfiguration();

    void initializeDataAnalysis();
    void doDataAnalysis();
};

#endif //__IPCPOSTPROCESS_HEADER_INCLUDED__
