#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "IPCsimulation.hpp"


void IPCsimulation::computeMSD() {
    double meanSquaredDisplacement = 0.0;
    for (IPC ipc: particles) {
        for (int j: {0, 1, 2}) {
            double delta_xj = ipc.ipcCenter.x[j] - ipcCentersPreviousPositions[ipc.number][j];
            relativePBC(delta_xj);
            displacementOfEachIPCs[ipc.number][j] += delta_xj;
            meanSquaredDisplacement += std::pow(displacementOfEachIPCs[ipc.number][j],2);
        }
    }
    meanSquaredDisplacement /= nIPCs;
    meanSquaredDisplFile << simulationTime*simulationTimeStep << "\t" << meanSquaredDisplacement << "\n";
}

void IPCsimulation::computeStaticProperties() {
    const std::vector<std::list<int>> listOfNeighbours = computeListOfBondedNeighbours();
    computeHistogramOfBondedNeighbours(listOfNeighbours);
}

std::vector<std::list<int>> IPCsimulation::computeListOfBondedNeighbours() {
    std::vector<std::list<int>> listOfNeighbours(nIPCs);

    for(int m=0; m<cells.getNumberofCells(); ++m)  // loop over all cells
    {
        const std::list<int> & ipcsInCurrentCell = cells.getIPCsInCell(m);
        const std::list<int> ipcsInNeighbouringCells = cells.getIPCsInNeighbouringCells(m);
        for(auto ipcInCell = ipcsInCurrentCell.cbegin(); ipcInCell != ipcsInCurrentCell.cend(); ++ipcInCell) {
            for(std::list<int>::const_iterator ipcInTheSameCell = std::next(ipcInCell); ipcInTheSameCell != ipcsInCurrentCell.cend(); ++ipcInTheSameCell) {
                double uij = computePotentialBetweenTwoIPCs(*ipcInCell, *ipcInTheSameCell);
                if (uij < 0.0) {
                    listOfNeighbours[*ipcInCell].push_back(*ipcInTheSameCell);
                    listOfNeighbours[*ipcInTheSameCell].push_back(*ipcInCell);
                }
            }
            for( auto ipcInTheOtherCells = ipcsInNeighbouringCells.cbegin(); ipcInTheOtherCells != ipcsInNeighbouringCells.cend(); ++ipcInTheOtherCells) {
                double uij = computePotentialBetweenTwoIPCs(*ipcInCell, *ipcInTheOtherCells);
                if (uij < 0.0) {
                    listOfNeighbours[*ipcInCell].push_back(*ipcInTheOtherCells);
                    listOfNeighbours[*ipcInTheOtherCells].push_back(*ipcInCell);
                }
            }
        }
    }

    /*
    for (int i = 0; i < nIPCs; ++i) {
        for (int j = i + 1; j < nIPCs; ++j) {
            double distance = 0.;
            for (int d: {0, 1, 2}) {
                double distance_d = particles[i].ipcCenter.x[d] - particles[j].ipcCenter.x[d];
                relativePBC(distance_d);
                distance += std::pow(distance_d, 2);
            }
            // if they are too far, it does not count
            if (distance > squaredInteractionRange)
                continue;
            //else
            localListOfNeighbours[i].push_back(j);
            localListOfNeighbours[j].push_back(i);
        }
    }*/
    return listOfNeighbours;
}

void IPCsimulation::computeHistogramOfBondedNeighbours(std::vector<std::list<int>> const& listOfNeighbours) {
    // compute how many bonded neighbours each particle has
    std::vector<int> numberOfNeighbours(nIPCs, 0);
    for (int i = 0; i < nIPCs; ++i)
        numberOfNeighbours[i] = listOfNeighbours[i].size();
    // compute the histogram of neighbours
    const int maxNumberOfNeighbours = *std::max_element(numberOfNeighbours.cbegin(), numberOfNeighbours.cend());
    std::vector<int> histogramOfNeighbours(maxNumberOfNeighbours+1, 0);
    for(int neighboursOfThisParticle: numberOfNeighbours)
        ++histogramOfNeighbours[neighboursOfThisParticle];

    // print and add to the averaged final histogram
    int i = 0;
    for(int n: histogramOfNeighbours) {
        numberOfNeighboursFile << simulationTime*simulationTimeStep << "\t" << i << "\t" << n << "\n";
        histogramOfBondedNeighbours[i] += n;

        ++i;
    }
}

double IPCsimulation::computePotentialBetweenTwoIPCs(const int firstIPC, const int secndIPC) {
    IPC const& first = particles[firstIPC];
    IPC const& secnd = particles[secndIPC];

    double binaryMixtureSign = 1.;
    if(binaryMixtureComposition > 0) {
        if ( (first.type == 'M'  && secnd.type == 'C') ||
             (first.type == 'C'  && secnd.type == 'M') ) {
            binaryMixtureSign = -1.;
        }
        else if( (first.type == 'M'  && secnd.type == 'M') ||
                   (first.type == 'C'  && secnd.type == 'C') ) {
        }
        else {
            std::cerr << __func__ << ": something really shitty is going on in the binary system recognition.";
            exit(1);
        }
    }

    // compute center-center distance
    double centerCenterSeparation[3];
    for (int i: {0, 1, 2}) {
        centerCenterSeparation[i] = first.ipcCenter.x[i] - secnd.ipcCenter.x[i];
        relativePBC(centerCenterSeparation[i]);
    }
    double centerCenterSeparationModulus = std::pow(centerCenterSeparation[0], 2)
                                         + std::pow(centerCenterSeparation[1], 2)
                                         + std::pow(centerCenterSeparation[2], 2);

    // if the CENTERS are too far, no interactions, skip this couple of IPCs
    if (centerCenterSeparationModulus >= squaredInteractionRange)
        return 0.0;

    // we are inside the interaction range; compute the interaction between centers
    double uij = 0.0;
    centerCenterSeparationModulus = std::sqrt(centerCenterSeparationModulus);
    const size_t centerCenterDistance = size_t( centerCenterSeparationModulus/forceAndEnergySamplingStep );
    uij += binaryMixtureSign*uBB[centerCenterDistance] + uHS[centerCenterDistance];

    // compute all the other 8 site-site separations
    double siteSiteSeparation[8][3];
    for (int i: {0, 1, 2}) {
        siteSiteSeparation[0][i] = first.ipcCenter.x[i] - secnd.firstPatch.x[i];
        siteSiteSeparation[1][i] = first.ipcCenter.x[i] - secnd.secndPatch.x[i];
        siteSiteSeparation[2][i] = first.firstPatch.x[i] - secnd.ipcCenter.x[i];
        siteSiteSeparation[3][i] = first.firstPatch.x[i] - secnd.firstPatch.x[i];
        siteSiteSeparation[4][i] = first.firstPatch.x[i] - secnd.secndPatch.x[i];
        siteSiteSeparation[5][i] = first.secndPatch.x[i] - secnd.ipcCenter.x[i];
        siteSiteSeparation[6][i] = first.secndPatch.x[i] - secnd.firstPatch.x[i];
        siteSiteSeparation[7][i] = first.secndPatch.x[i] - secnd.secndPatch.x[i];
        for (int j = 0; j < 8; ++j)
            relativePBC(siteSiteSeparation[j][i]);
    }

    // compute all the other 8 site-site interactions
    for (int j = 0; j < 8; ++j) {
        double siteSiteSeparationModulus = siteSiteSeparation[j][0]*siteSiteSeparation[j][0]
                                         + siteSiteSeparation[j][1]*siteSiteSeparation[j][1]
                                         + siteSiteSeparation[j][2]*siteSiteSeparation[j][2];

        // if we are too far, no interaction, skip to the next site-site pair
        if (siteSiteSeparationModulus >= squaredInteractionRange)
            continue;

        siteSiteSeparationModulus = std::sqrt(siteSiteSeparationModulus);
        const size_t dist = size_t( siteSiteSeparationModulus/forceAndEnergySamplingStep );
        if (j == 0) { // center - patch1
            uij += binaryMixtureSign*uBs1[dist];
        } else if (j == 1) { // center - patch2
            uij += binaryMixtureSign*uBs2[dist];
        } else if (j == 2) { // patch1 - center
            uij += binaryMixtureSign*uBs1[dist];
        } else if (j == 3) { // patch1 - patch1
            uij += binaryMixtureSign*us1s1[dist];
        } else if (j == 4) { // patch1 - patch2
            uij += binaryMixtureSign*us1s2[dist];
        } else if (j == 5) { // patch2 - center
            uij += binaryMixtureSign*uBs2[dist];
        } else if (j == 6) { // patch2 - patch1
            uij += binaryMixtureSign*us1s2[dist];
        } else if (j == 7) { // patch2 - patch2
            uij += binaryMixtureSign*us2s2[dist];
        } else {
            std::cerr << __func__ << ":: something really shitty is going on in the pair potential computation.";
            exit(1);
        }
    }
    return uij;
}

void IPCsimulation::printFinalHistogramOfBondedNeighbours() {
    double norm = printingInterval/simulationTotalDuration;
    std::ofstream averageNumberOfNeighboursFile("siml/averageNumberOfNeighbours.out");
    averageNumberOfNeighboursFile << std::scientific << std::setprecision(6);

    int i = 0;
    for(int n: histogramOfBondedNeighbours) {
        averageNumberOfNeighboursFile << i << "\t" << n*norm << "\n";
        ++i;
    }
}