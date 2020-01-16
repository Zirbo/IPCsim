#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <unordered_set>
#include "IPCsimulation.hpp"



void IPCsimulation::initializeDataAnalysis() {
    pairCorrelation.initialize(20, simulationBoxSide, nIPCs);

    meanSquaredDisplFile.open("siml/meanSquaredDisplacement.out");
    meanSquaredDisplFile << std::scientific << std::setprecision(6);
    ipcCentersPreviousPositions.resize(nIPCs, {0.0, 0.0, 0.0});
    displacementOfEachIPCs.resize(nIPCs, {0.0, 0.0, 0.0});

    numberOfNeighboursFile.open("siml/numberOfNeighbours.out");
    numberOfNeighboursFile << std::scientific << std::setprecision(6);

    clusterSizesFile.open("siml/clusterSizes.out");
    clusterSizesFile << std::scientific << std::setprecision(6);

    ipcOrientations.resize(nIPCs, {0.0, 0.0, 0.0});
    nematicOrderParameter.resize(nIPCs, 0.0);
    nematicOrderParameterFile.open("siml/NOP.out");
    nematicOrderParameterFile << std::scientific << std::setprecision(6);
}

void IPCsimulation::doDataAnalysis() {
    const std::vector<std::list<int>> listOfNeighbours = computeListOfBondedNeighbours();
    computeHistogramOfBondedNeighbours(listOfNeighbours);
    computeClusters(listOfNeighbours);

    updateOrientations();
    computeNematicOrderParameter(listOfNeighbours);
}

void IPCsimulation::printDataAnalysis() {
    const double g_r_integral = pairCorrelation.print("siml/g_r.out");
    outputFile << "The integral of g(r) is " << g_r_integral << " and is should be equal to the number of particles minus one, " << nIPCs-1 << std::endl;
    printHistogramOfBondedNeighbours();
    printClusterSizes();
    printNematicOrderPatameter();
}

/******************************************************************************
 ****************DYNAMIC PROPERTIES*******************************************
 ******************************************************************************/

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

void IPCsimulation::updateOrientations() {
    const double norm = 1./firstPatchEccentricity;
    for (IPC ipc: particles) {
        for (int d: {0, 1, 2}) {
            ipcOrientations[ipc.number][d] = ipc.firstPatch.x[d] - ipc.ipcCenter.x[d];
            relativePBC(ipcOrientations[ipc.number][d]);
            ipcOrientations[ipc.number][d] *= norm;
        }
    }
}
void IPCsimulation::initializeAutocorrelations() {

}
void IPCsimulation::computeAutocorrelations() {

}

/******************************************************************************
 ****************STATIC PROPERTIES*********************************************
 ******************************************************************************/

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

    std::map<int, std::unordered_set<int>> clusters;

    for (int i = 0; i < nIPCs; ++i) {
        std::unordered_set<int> currentGroup;
        currentGroup.insert(i);
        for (int j: listOfNeighbours[i]) {
            currentGroup.insert(j);
        }
        // find if any IPC of this group is present in any enstablished cluster
        std::unordered_set<int> matches;

        for (int ipc: currentGroup) {
            for (auto cluster: clusters) {
                if (cluster.second.count(ipc) == 1) {
                    // match found!
                    matches.insert(cluster.first);
                }
            }
        }
        // if any match has been found, add the current group to the smallest head of chain;
        // if more than one match is present, merge to the smallest head of chain all the others
        if(!matches.empty()) {
            // find the smallest head of chain
            int minimum = *std::min_element(matches.cbegin(), matches.cend());
            // added the current group
            for (int ipc: currentGroup)
                clusters[minimum].insert(ipc);
            // add all the others
            for (auto head: matches) {
                if (head == minimum)
                    continue;
                for (int ipc: clusters[head])
                    clusters[minimum].insert(ipc);
                clusters.erase(head);
            }
        }
        else {
            // create new cluster
            clusters.insert( std::make_pair(i, currentGroup) );
        }
    }


    // analyze the clusters
    std::map<int, int> localClusterSizes;
    for (auto i: clusters) {
        const int size = i.second.size();
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
        clusterSizesFile << simulationTime*simulationTimeStep << "\t" << n.first << "\t" << n.second << "\n";

        integral += n.first*n.second;
        if (clusterSizes.count(n.first) == 0) {
            clusterSizes[n.first] = n.second;
        }
        else {
            clusterSizes[n.first] += n.second;
        }
    }
    if(integral != nIPCs) {
        std::cerr << __func__ << ": something really shitty is going on in the cluster size analysis.";
        exit(1);
    }
}

void IPCsimulation::printClusterSizes() {
    double norm = printingInterval/simulationTotalDuration;
    std::ofstream clusterSizesFile("siml/averageClusterSizes.out");
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

void IPCsimulation::computeNematicOrderParameter(std::vector<std::list<int>> const& listOfNeighbours) {
    // loop on the lists of neighbours
    double averageNOP = 0.0;
    for (int i = 0; i < nIPCs; ++i) {
        double modulusNOPi = 0.;
        // loop on the neighbours inside the list
        for (int j: listOfNeighbours[i]) {
            double modulusNOPij = 0.;
            for (int d: {0, 1, 2}) {
                modulusNOPij += ipcOrientations[i][d]*ipcOrientations[j][d];
            }
            modulusNOPi += std::pow(modulusNOPij,2);
        }
        // normalize and print
        if(!listOfNeighbours[i].empty())
            modulusNOPi = 1.5*modulusNOPi/listOfNeighbours[i].size() - 0.5;
        nematicOrderParameterFile << simulationTime*simulationTimeStep << "\t" << i << "\t" << modulusNOPi << "\n";
        nematicOrderParameter[i] += modulusNOPi;
        averageNOP += modulusNOPi;
    }
    nematicOrderParameterFile << "# Global average of the nematic order parameter at time "
                              << simulationTime*simulationTimeStep <<": " << averageNOP/nIPCs << "!\n";
}

void IPCsimulation::printNematicOrderPatameter() {
    double norm = printingInterval/simulationTotalDuration;
    std::ofstream nematicOrderParameterFile("siml/averageNOP.out");
    nematicOrderParameterFile << std::scientific << std::setprecision(6);

    double averageNOP = 0.0;
    for (int i = 0; i < nIPCs; ++i) {
        nematicOrderParameterFile << i << "\t" << nematicOrderParameter[i]*norm << "\n";
        averageNOP += nematicOrderParameter[i];
    }
    averageNOP *= norm;
    nematicOrderParameterFile << "# Average of the final NOP: " << averageNOP/nIPCs << std::endl;
    nematicOrderParameterFile.close();
}
