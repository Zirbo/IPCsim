#include <iostream>
#include <fstream>
#include "IPCpostprocessNeighbourAnalysis.hpp"

void IPCneighboursAnalysis::accumulate(IPCpotential const& potential, Ensemble const & ipcs) {
    auto listOfBondedNeighbours = computeListOfBondedNeighbours(potential, ipcs);
    computeHistogramOfBondedNeighbours(listOfBondedNeighbours, ipcs);
}

std::vector<std::list<int>> IPCneighboursAnalysis::computeListOfBondedNeighbours(IPCpotential const& potential, Ensemble const & ipcs) {
    const int nIPCs = (int)ipcs.size();
    std::vector<std::list<int>> listOfBondedNeighbours(nIPCs);

    for (int i = 0; i < nIPCs; ++i) {
        for (int j = i + 1; j < nIPCs - 1; ++j) {
            //
            double uij = computePotentialBetweenTwoIPCs(potential, ipcs[i], ipcs[j]);
            if (uij < 0.0) {
                listOfBondedNeighbours[i].push_back(j);
                listOfBondedNeighbours[j].push_back(i);
            }
        }
    }
    return listOfBondedNeighbours;
}

double IPCneighboursAnalysis::computePotentialBetweenTwoIPCs(IPCpotential const& potential, IPC const& firstIPC, IPC const& secndIPC) {

    // compute center-center distance
    double centerCenterDistance = 0.;
    for (int d: DIMENSIONS) {
        double centerCenterSeparation = firstIPC.ipcCenter.x[d] - secndIPC.ipcCenter.x[d];
        relativePBC(centerCenterSeparation);
        centerCenterDistance += std::pow(boxSide[d]*centerCenterSeparation, 2);
    }
    centerCenterDistance = std::sqrt(centerCenterDistance);

    // if the CENTERS are too far, no interactions, skip this couple of IPCs
    if (centerCenterDistance >= interactionRange)
        return 0.0;

    // we are inside the interaction range; compute the interaction between centers
    double uij = 0.0;
    const size_t dist = size_t( centerCenterDistance/potential.spacing );
    uij += potential.uBB[dist];

    // compute all the other 8 site-site separations
    double siteSiteSeparation[8][3];
    for (int d: DIMENSIONS) {
        siteSiteSeparation[0][d] = firstIPC.ipcCenter.x[d]  - secndIPC.firstPatch.x[d];
        siteSiteSeparation[1][d] = firstIPC.ipcCenter.x[d]  - secndIPC.secndPatch.x[d];
        siteSiteSeparation[2][d] = firstIPC.firstPatch.x[d] - secndIPC.ipcCenter.x[d];
        siteSiteSeparation[3][d] = firstIPC.firstPatch.x[d] - secndIPC.firstPatch.x[d];
        siteSiteSeparation[4][d] = firstIPC.firstPatch.x[d] - secndIPC.secndPatch.x[d];
        siteSiteSeparation[5][d] = firstIPC.secndPatch.x[d] - secndIPC.ipcCenter.x[d];
        siteSiteSeparation[6][d] = firstIPC.secndPatch.x[d] - secndIPC.firstPatch.x[d];
        siteSiteSeparation[7][d] = firstIPC.secndPatch.x[d] - secndIPC.secndPatch.x[d];
        for (int j = 0; j < 8; ++j)
            relativePBC(siteSiteSeparation[j][d]);
    }

    // compute all the other 8 site-site interactions
    for (int j = 0; j < 8; ++j) {
        double siteSiteDistance = 0.;
        for (int d: DIMENSIONS)
            siteSiteDistance += std::pow(boxSide[d]*siteSiteSeparation[j][d], 2);
        siteSiteDistance = std::sqrt(siteSiteDistance);

        // if we are too far, no interaction, skip to the next site-site pair
        if (siteSiteDistance >= interactionRange)
            continue;

        const size_t dist = size_t( siteSiteDistance/potential.spacing );
        if (j == 0 || j == 1 || j == 2 || j == 5) {
            uij += potential.uBs[dist];
        } else if (j == 3 || j == 4 || j == 6 || j == 7) {
            uij += potential.uss[dist];
        } else {
            std::cerr << __func__ << ":: something bad happened!";
            exit(1);
        }
    }
    return uij;
}

void IPCneighboursAnalysis::computeHistogramOfBondedNeighbours(std::vector<std::list<int>> const& listOfNeighbours, Ensemble const& ipcs) {
    // compute how many bonded neighbours each particle has
    int nIPCs = ipcs.size();
    std::vector<int> numberOfNeighbours(nIPCs, 0);
    for (int i = 0; i < nIPCs; ++i)
        numberOfNeighbours[i] = listOfNeighbours[i].size();

    // compute the histogram of neighbours
    std::map<int, int> histogramOfNeighbours;
    for(int neighboursOfThisParticle: numberOfNeighbours) {
        if (histogramOfNeighbours.count(neighboursOfThisParticle) == 0) {
            histogramOfNeighbours[neighboursOfThisParticle] = 1;
        }
        else {
            histogramOfNeighbours[neighboursOfThisParticle] += 1;
        }
    }
    // print and add to the averaged final histogram
    for(auto n: histogramOfNeighbours) {

        if (histogramOfBondedNeighbours.count(n.first) == 0) {
            histogramOfBondedNeighbours[n.first] = n.second;
        }
        else {
            histogramOfBondedNeighbours[n.first] += n.second;
        }
    }
    totalSamples++;
}

void IPCneighboursAnalysis::print(std::string const& outputFileName) {
    const double norm = 1./totalSamples;
    std::ofstream averageNumberOfNeighboursFile(outputFileName);
    //averageNumberOfNeighboursFile << std::scientific << std::setprecision(2);

    for(auto n: histogramOfBondedNeighbours) {
        averageNumberOfNeighboursFile << n.first << "\t" << norm*n.second << "\n";
    }
}
