#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <set>
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
    computeClusters(listOfNeighbours);
}

void IPCsimulation::printStaticProperties() {
    const double g_r_integral = pairCorrelation.print("siml/g_r.out");
    outputFile << "The integral of g(r) is " << g_r_integral << " and is should be equal to the number of particles minus one, " << nIPCs-1 << std::endl;
    printHistogramOfBondedNeighbours();
    printClusterSizes();
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
    return listOfNeighbours;
}

void IPCsimulation::computeHistogramOfBondedNeighbours(std::vector<std::list<int>> const& listOfNeighbours) {
    // compute how many bonded neighbours each particle has
    std::vector<int> numberOfNeighbours(nIPCs, 0);
    for (int i = 0; i < nIPCs; ++i)
        numberOfNeighbours[i] = listOfNeighbours[i].size();
    // compute the histogram of neighbours
    std::map<int, int> histogramOfNeighbours;
    //const int maxNumberOfNeighbours = *std::max_element(numberOfNeighbours.cbegin(), numberOfNeighbours.cend());
    //std::vector<int> histogramOfNeighbours(maxNumberOfNeighbours+1, 0);

    for(int neighboursOfThisParticle: numberOfNeighbours) {
        if (histogramOfNeighbours.count(neighboursOfThisParticle) == 0) {
            histogramOfNeighbours[neighboursOfThisParticle] = 1;
        }
        else {
            histogramOfNeighbours[neighboursOfThisParticle] += 1;
        }
    }

        //++histogramOfNeighbours[neighboursOfThisParticle];

    // print and add to the averaged final histogram
    for(auto n: histogramOfNeighbours) {
        numberOfNeighboursFile << simulationTime*simulationTimeStep << "\t" << n.first << "\t" << n.second << "\n";

        if (histogramOfBondedNeighbours.count(n.first) == 0) {
            histogramOfBondedNeighbours[n.first] = n.second;
        }
        else {
            histogramOfBondedNeighbours[n.first] += n.second;
        }
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

void IPCsimulation::printHistogramOfBondedNeighbours() {
    double norm = printingInterval/simulationTotalDuration;
    std::ofstream averageNumberOfNeighboursFile("siml/averageNumberOfNeighbours.out");
    averageNumberOfNeighboursFile << std::scientific << std::setprecision(6);

    for(auto n: histogramOfBondedNeighbours) {
        averageNumberOfNeighboursFile << n.first << "\t" << norm*n.second << "\n";
    }
}

void IPCsimulation::computeClusters(std::vector<std::list<int>> const& listOfNeighbours) {
    // /home/bianchi/IPC-QUASI-2D-MC/IPC-Quasi2D-PostProcessing-NEW/ipc-postprocessing.f90

    // subdivide the IPCs into clusters
    std::list<std::set<int>> clusters;
    for (int i = 0; i < nIPCs; ++i) {
        bool aClusterWasFound = false;
        std::set<int>* clusterMatch;

        // is this IPC already in a cluster? if yes, skip it
        for (auto &cl: clusters) {
            if (cl.count(i) == 1) {
                aClusterWasFound = true;
                clusterMatch = &cl;
                break;
            }
        }

        // check if any of its bonded neighbours belong to an existing cluster
        for (auto j: listOfNeighbours[i]) {
            // find in current clusters
            for (auto &cl: clusters) {
                if (cl.count(j) == 1) {
                    if (!aClusterWasFound) {
                        // found a first match with a cluster
                        cl.insert(i); // this line is fucking up everything...
                        aClusterWasFound = true;
                        clusterMatch = &cl;
                    }
                    else {
                        // found a new match; we need to merge the cluster of the first match with the one we just found
                        clusterMatch->insert(cl.begin(), cl.end());
                        // then clean the new one
                        cl.clear();
                    }
                }
            }
        }
        if (!aClusterWasFound) {
            // no cluster was found at all; let's create a new one
            std::set<int> a;
            a.insert(i);
            for (auto j: listOfNeighbours[i])
                a.insert(j);
            clusters.push_back(a);
        }
    }

    clusters.remove(std::set<int>());


    // analyze the clusters
    std::map<int, int> localClusterSizes;
    for (auto i: clusters) {
        const int size = i.size();
        if (localClusterSizes.count(size) == 0) {
            localClusterSizes[size] = 1;
        }
        else {
            localClusterSizes[size] += 1;
        }
    }

    // print and add to the averaged final histogram
    int integral = 0;
    for(std::pair<int, int> n: localClusterSizes) {
        integral += n.first*n.second;
        if (clusterSizes.count(n.first) == 0) {
            clusterSizes[n.first] = n.second;
        }
        else {
            clusterSizes[n.first] += n.second;
        }
    }
    if(integral != nIPCs) {
        std::cerr << "error";
        exit(1);
    }
}

void IPCsimulation::printClusterSizes() {
    double norm = printingInterval/simulationTotalDuration;
    std::ofstream clusterSizesFile("siml/clusterSizes.out");
    clusterSizesFile << std::scientific << std::setprecision(6);

    int i = 0;
    int integral = 0;
    for(auto n: clusterSizes) {
        clusterSizesFile << n.first << "\t" << norm*n.second << "\n";
        integral += n.first*n.second;
        ++i;
    }
    clusterSizesFile << "#" << integral << " = " << nIPCs;
    clusterSizesFile.close();
}
