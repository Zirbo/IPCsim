#ifndef __IPCPOSTPROCESS_DATAANALYSIS_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_DATAANALYSIS_HEADER_INCLUDED__

#include <map>
#include <vector>
#include <list>
#include <cmath>

#include "IPC.hpp"
#include "IPCpostprocessPotential.hpp"

class IPCneighboursAnalysis {
public:
    IPCneighboursAnalysis(double pBoxSideX, double pBoxSideY, double pBoxSideZ, double pInteractionRange)
        : boxSideX(pBoxSideX), boxSideY(pBoxSideY), boxSideZ(pBoxSideZ),
          interactionRange(pInteractionRange), totalSamples(0) {}
    void accumulate(IPCpotential const& potential, std::vector<IPC> const & particles);
    void print();

private:
    std::map<int, int> histogramOfBondedNeighbours;
    //std::ofstream numberOfNeighboursFile;
    double boxSideX, boxSideY, boxSideZ;
    double interactionRange;
    int totalSamples;

    std::vector<std::list<int>> computeListOfBondedNeighbours(IPCpotential const& potential, std::vector<IPC> const& particles);
    inline void relativePBC(double & x) {  x -= std::round(x);  }
    double computePotentialBetweenTwoIPCs(IPCpotential const& pPotential, IPC const& firstIPC, IPC const& secndIPC);
    void computeHistogramOfBondedNeighbours(std::vector<std::list<int>> const& listOfNeighbours, std::vector<IPC> const& particles);
    void printHistogramOfBondedNeighbours();

};

#endif //__IPCPOSTPROCESS_DATAANALYSIS_HEADER_INCLUDED__
