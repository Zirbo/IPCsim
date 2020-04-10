#include <cmath>
#include <fstream>
#include <iomanip>

#include "IPCpostprocessOrientationsAnalysis.hpp"

void IPCorientationsAnalysis::accumulate(SpaceVector const& ipcOrientations) {
    accumulateOrientationsHistogram(ipcOrientations);
}
void IPCorientationsAnalysis::print() {
    printOrientationsHistogram();
}

void IPCorientationsAnalysis::accumulateOrientationsHistogram(SpaceVector const& ipcOrientations) {
    // polar -> angle with the z-axis; azimuth -> angle with the x-axis
    // accumulute in a histogram the polar and azimuth angles of each IPC
    // we assume patch symmetry, so at the moment if cosPolar < 0 we do polar += pi/2 and azimuth += pi
    const double binSize = M_PI/orientationHistogramSize;

    for (auto const& ipcOrientation: ipcOrientations) {
        double polarAngle = std::acos(ipcOrientation[2]);
        const double polarScaling = 1./std::sin(polarAngle);
        double azimuthAngle = M_PI + std::atan2(ipcOrientation[1]*polarScaling, ipcOrientation[0]*polarScaling);

/*        if (ipcOrientation[2] < 0) {
            polarAngle = M_PI - polarAngle;
            azimuthAngle += M_PI;
            if (azimuthAngle > 2*M_PI)
                azimuthAngle -= 2*M_PI;
        }*/

        const int polarBin = int(polarAngle/binSize);
        const int azimuthBin = int(azimuthAngle/binSize);
        ++orientationsHistogram[polarBin][azimuthBin];
    }

    ++totalSamples;
}

void IPCorientationsAnalysis::printOrientationsHistogram() {
    std::ofstream outputFile("orientationsHistogram.out");
    outputFile << std::scientific << std::setprecision(6);

    double norm = 1./totalSamples;
    const double binSize = M_PI/orientationHistogramSize;

    for (int a = 0; a < 2*orientationHistogramSize; ++a) {
        for (int p = 0; p < orientationHistogramSize; ++p) {
            outputFile << a*binSize << "\t" << p*binSize << "\t" << orientationsHistogram[p][a]*norm << "\t" << std::log10(orientationsHistogram[p][a]*norm) << "\n";
        }
    }
    outputFile.close();
}
