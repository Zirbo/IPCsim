#include <iostream>
#include <fstream>
#include "IPCpostprocessNeighbourAnalysis.hpp"

void IPCneighboursAnalysis::accumulate(IPCpotential const& potential, Ensemble const & ipcs) {
    auto listOfBondedNeighbours = computeListOfBondedNeighbours(potential, ipcs);
    computeHistogramOfBondedNeighbours(listOfBondedNeighbours, ipcs);
}

void IPCneighboursAnalysis::print() {
    printHistogramOfBondedNeighbours();
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
    double centerCenterSeparation[3];
    for (int i: {0, 1, 2}) {
        centerCenterSeparation[i] = firstIPC.ipcCenter.x[i] - secndIPC.ipcCenter.x[i];
        relativePBC(centerCenterSeparation[i]);
    }
    double centerCenterSeparationModulus = std::pow(boxSideX*centerCenterSeparation[0], 2)
                                         + std::pow(boxSideY*centerCenterSeparation[1], 2)
                                         + std::pow(boxSideZ*centerCenterSeparation[2], 2);
    centerCenterSeparationModulus = std::sqrt(centerCenterSeparationModulus);

    // if the CENTERS are too far, no interactions, skip this couple of IPCs
    if (centerCenterSeparationModulus >= interactionRange)
        return 0.0;

    // we are inside the interaction range; compute the interaction between centers
    double uij = 0.0;
    const size_t centerCenterDistance = size_t( centerCenterSeparationModulus/potential.spacing );
    uij += potential.uBB[centerCenterDistance];

    // compute all the other 8 site-site separations
    double siteSiteSeparation[8][3];
    for (int i: {0, 1, 2}) {
        siteSiteSeparation[0][i] = firstIPC.ipcCenter.x[i]  - secndIPC.firstPatch.x[i];
        siteSiteSeparation[1][i] = firstIPC.ipcCenter.x[i]  - secndIPC.secndPatch.x[i];
        siteSiteSeparation[2][i] = firstIPC.firstPatch.x[i] - secndIPC.ipcCenter.x[i];
        siteSiteSeparation[3][i] = firstIPC.firstPatch.x[i] - secndIPC.firstPatch.x[i];
        siteSiteSeparation[4][i] = firstIPC.firstPatch.x[i] - secndIPC.secndPatch.x[i];
        siteSiteSeparation[5][i] = firstIPC.secndPatch.x[i] - secndIPC.ipcCenter.x[i];
        siteSiteSeparation[6][i] = firstIPC.secndPatch.x[i] - secndIPC.firstPatch.x[i];
        siteSiteSeparation[7][i] = firstIPC.secndPatch.x[i] - secndIPC.secndPatch.x[i];
        for (int j = 0; j < 8; ++j)
            relativePBC(siteSiteSeparation[j][i]);
    }

    // compute all the other 8 site-site interactions
    for (int j = 0; j < 8; ++j) {
        double siteSiteSeparationModulus = std::pow(boxSideX*siteSiteSeparation[j][0], 2)
                                         + std::pow(boxSideY*siteSiteSeparation[j][1], 2)
                                         + std::pow(boxSideZ*siteSiteSeparation[j][2], 2);
        siteSiteSeparationModulus = std::sqrt(siteSiteSeparationModulus);

        // if we are too far, no interaction, skip to the next site-site pair
        if (siteSiteSeparationModulus >= interactionRange)
            continue;

        const size_t dist = size_t( siteSiteSeparationModulus/potential.spacing );
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

void IPCneighboursAnalysis::printHistogramOfBondedNeighbours() {
    const double norm = 1./totalSamples;
    std::ofstream averageNumberOfNeighboursFile("numberOfBondedNeighbours.out");
    //averageNumberOfNeighboursFile << std::scientific << std::setprecision(2);

    for(auto n: histogramOfBondedNeighbours) {
        averageNumberOfNeighboursFile << n.first << "\t" << norm*n.second << "\n";
    }
}
