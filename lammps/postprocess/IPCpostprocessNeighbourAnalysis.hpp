#ifndef __IPCPOSTPROCESS_NEIGHBOURSANALYSIS_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_NEIGHBOURSANALYSIS_HEADER_INCLUDED__

#include <map>
#include <vector>
#include <list>
#include <cmath>

#include "IPC.hpp"
#include "IPCpostprocessPotential.hpp"

class IPCneighboursAnalysis {
public:
    IPCneighboursAnalysis(Triad pBoxSides, double pInteractionRange)
        : boxSide{pBoxSides},
          interactionRange(pInteractionRange), totalSamples(0) {}
    void accumulate(IPCpotential const& potential, const Ensemble &ipcs);
    void print(std::string const& outputFileName);

private:
    IPCneighboursAnalysis();
    std::map<int, int> histogramOfBondedNeighbours;
    Triad boxSide;
    double interactionRange;
    int totalSamples;

    std::vector<std::list<int>> computeListOfBondedNeighbours(IPCpotential const& potential, const Ensemble &ipcs);
    inline void relativePBC(double & x) {  x -= std::lround(x);  }
    double computePotentialBetweenTwoIPCs(IPCpotential const& pPotential, IPC const& firstIPC, IPC const& secndIPC);
    void computeHistogramOfBondedNeighbours(std::vector<std::list<int>> const& listOfNeighbours, const Ensemble &ipcs);

};

#endif //__IPCPOSTPROCESS_NEIGHBOURSANALYSIS_HEADER_INCLUDED__
