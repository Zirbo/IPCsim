#ifndef __IPCPOSTPROCESS_ORIENTATIONSANALYSIS_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_ORIENTATIONSANALYSIS_HEADER_INCLUDED__

#include <vector>
//#include <cmath>

#include "IPC.hpp"
#include "IPCpostprocessPotential.hpp"

class IPCorientationsAnalysis {
public:
    IPCorientationsAnalysis()
        : orientationHistogramSize(40) {
        orientationsHistogram.resize(orientationHistogramSize, std::vector<double>(2*orientationHistogramSize, 0.));
    }
    void accumulate(SpaceVector const& ipcOrientations);
    void print();

private:
    std::vector<std::vector<double>> orientationsHistogram;
    int orientationHistogramSize;
    int totalSamples;

    inline void relativePBC(double & x) {  x -= std::round(x);  }

    void accumulateOrientationsHistogram(SpaceVector const& ipcOrientations);
    void printOrientationsHistogram();

};

#endif //__IPCPOSTPROCESS_ORIENTATIONSANALYSIS_HEADER_INCLUDED__
