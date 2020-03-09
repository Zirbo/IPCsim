#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "IPCsimulation.hpp"



void IPCsimulation::initializeDataAnalysis() {
    pairCorrelation.initialize(40, simulationBoxSide, nIPCs);

    meanSquaredDisplFile.open("siml/meanSquaredDisplacement.out");
    meanSquaredDisplFile << std::scientific << std::setprecision(6);
    ipcCentersPreviousPositions.resize(nIPCs, {0.0, 0.0, 0.0});
    updatePreviousPositions();
    displacementOfEachIPCs.resize(nIPCs, {0.0, 0.0, 0.0});
    computeMSD();

    numberOfNeighboursFile.open("siml/numberOfNeighbours.out");
    numberOfNeighboursFile << std::scientific << std::setprecision(6);

    clusterSizesFile.open("siml/clusterSizes.out");
    clusterSizesFile << std::scientific << std::setprecision(6);
    average_pOverL = 0.;

    initializeAutocorrelations();

    nematicOrderParameter.resize(nIPCs, 0.0);
    nematicOrderParameterFile.open("siml/NOP.out");
    nematicOrderParameterFile << std::scientific << std::setprecision(6);

    orientationHistogramSize = 40;
    orientationsHistogram.resize(orientationHistogramSize, std::vector<double>(2*orientationHistogramSize, 0.));
}

void IPCsimulation::doDataAnalysis() {
    computeMSD();
    updateOrientations();
    computeAutocorrelations();

    pairCorrelation.compute(particles);
    const std::vector<std::list<int>> listOfNeighbours = computeListOfBondedNeighbours();
    computeHistogramOfBondedNeighbours(listOfNeighbours);
    doClustersAnalysis(listOfNeighbours);
    computeNematicOrderParameter(listOfNeighbours);

    accumulateOrientationsHistogram();

    updatePreviousPositions();
}

void IPCsimulation::printDataAnalysis() {
    const double g_r_integral = pairCorrelation.print("siml/g_r.out");
    outputFile << "The integral of g(r) is " << g_r_integral << " and it should be equal to the number of particles minus one, " << nIPCs-1 << std::endl;
    printHistogramOfBondedNeighbours();
    printClusterSizes();
    printNematicOrderPatameter();
    printOrientationsHistogram();
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

void IPCsimulation::updatePreviousPositions() {
    // set the current as the previous iteration state
    for (IPC ipc: particles) {
        for (int j: {0, 1, 2}) {
            ipcCentersPreviousPositions[ipc.number][j] = ipc.ipcCenter.x[j];
        }
    }
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
    ipcOrientations.resize(nIPCs, {0.0, 0.0, 0.0});
    initialOrientations.resize(nIPCs, {0.0, 0.0, 0.0});
    initialVelocities.resize(nIPCs, {0.0, 0.0, 0.0});
    updateOrientations();

    normOfCv = 0.0;
    for(int i = 0; i < nIPCs; ++i) {
        for (int d: {0, 1, 2}) {
            initialOrientations[i][d] = ipcOrientations[i][d];
            initialVelocities[i][d] = particles[i].ipcCenter.v[d];
            normOfCv += std::pow(initialVelocities[i][d],2);
        }
    }
    normOfCv = 1./normOfCv;
    normOfCn = 1./nIPCs;

    autocorrelationsFile.open("siml/autocorrelations.out");
    autocorrelationsFile << std::scientific << std::setprecision(6);

    computeAutocorrelations();
}

void IPCsimulation::computeAutocorrelations() {
    double Cn = 0.0;
    double Cv = 0.0;
    for(int i = 0; i < nIPCs; ++i) {
        for (int d: {0, 1, 2}) {
            Cn += ipcOrientations[i][d]*initialOrientations[i][d];
            Cv += particles[i].ipcCenter.v[d]*initialVelocities[i][d];
        }
    }
    Cn *= normOfCn;
    Cv *= normOfCv;

    autocorrelationsFile << simulationTime*simulationTimeStep << "\t"
                         << Cn << "\t" << Cv << "\n";
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

    // if requested, override the IPC type with the NoN so that it gets printed
    if(overrideTypeWithNumberOfNeighbours) {
        for(int i = 0; i < nIPCs; ++i) {
            particles[i].type = 'A' + numberOfNeighbours[i];
        }
    }

    // if requested, override the IPC type with the NoN so that it gets printed
    if(overrideTypeWithNumberOfNeighboursThreshold != 0) {
        for(int i = 0; i < nIPCs; ++i) {
            char newType = (numberOfNeighbours[i] > overrideTypeWithNumberOfNeighboursThreshold)? 'W' : 'S';
            particles[i].type = newType;
        }
    }

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

void IPCsimulation::doClustersAnalysis(std::vector<std::list<int>> const& listOfNeighbours) {
    // find and count the clusters
    std::map<int, std::unordered_set<int>> clusters = findClusters(listOfNeighbours);
    std::map<int, int> localClusterSizesHistogram = computeClusterSizesHistogram(clusters);

    // print and add to the cluster size histogram that averages on all simulation
    int integral = 0;
    for(std::pair<int, int> n: localClusterSizesHistogram) {
        clusterSizesFile << simulationTime*simulationTimeStep << "\t" << n.first << "\t" << n.second << "\n";

        integral += n.first*n.second;
        if (clusterSizesHistogram.count(n.first) == 0) {
            clusterSizesHistogram[n.first] = n.second;
        }
        else {
            clusterSizesHistogram[n.first] += n.second;
        }
    }
    if(integral != nIPCs) {
        std::cerr << __func__ << ": something really shitty is going on in the cluster size analysis.";
        exit(1);
    }

    // these are only for chains, no point in computing it for any other system -- need to improve the flag
    chainFlatnessAnalysis(listOfNeighbours, clusters);
    if(overrideTypeWithClusterID)
        overrideIPCtypeWithClusterID(clusters);
}

std::map<int, std::unordered_set<int>> IPCsimulation::findClusters(std::vector<std::list<int>> const& listOfNeighbours) {
    // first is the cluster head-of-chain IPC, second is the unordered set of the IPCs in the cluster, INCLUDING the head
    std::map<int, std::unordered_set<int>> clusters;

    /* Interate through the IPCs, and for each make a group with all its neighbours.
     * Check if any of them is present in other clusters; if so, join all the clusters, otherwise create a new one.
     * The join is done in two steps: first find all the matches, and only afterwards to the actual joining.
     * This is to avoid changing the map over which we are iterating.
     */
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

    return clusters;
}
std::map<int, int> IPCsimulation::computeClusterSizesHistogram(std::map<int, std::unordered_set<int>> const& clusters) {
    std::map<int, int> clustersSizeHistogram;
    for (auto i: clusters) {
        const int size = i.second.size();
        if (clustersSizeHistogram.count(size) == 0) {
            clustersSizeHistogram[size] = 1;
        }
        else {
            clustersSizeHistogram[size] += 1;
        }
    }

    return clustersSizeHistogram;
}

void IPCsimulation::chainFlatnessAnalysis(std::vector<std::list<int>> const& listOfNeighbours, std::map<int, std::unordered_set<int>> const& clusters) {
    // go through all the clusters, compute the end-to-end distance, and divide it by the cluster size. Average this number through all the simulation.
    double pOverL = 0.0;
    int pOverLsamples = 0;

    //const size_t maxSize = size_t(simulationBoxSide*std::cbrt(3.));

    for (auto cluster: clusters) {
        const size_t clusterSize = cluster.second.size();
        if ( clusterSize == 0 ) // || clusterSize < maxSize )
            continue;
        // find the endpoints of this cluster --- watch out for branchpoints :P
        int firstEndpoint = -1;
        int secondEndpoint = -1;
        for (int IPC: cluster.second) {
            if(listOfNeighbours[IPC].size() == 1) {
                if (firstEndpoint == -1)
                    firstEndpoint = IPC;
                else if (secondEndpoint == -1) {
                    secondEndpoint = IPC;
                    break;
                }
            }
        }
        // compute the distance between endpoints and divide for the cluster size!
        if (firstEndpoint != -1 && secondEndpoint != -1) {
            double endToEndDistance[3];
            for (int i: {0, 1, 2}) {
                endToEndDistance[i] = particles[firstEndpoint].ipcCenter.x[i] - particles[secondEndpoint].ipcCenter.x[i];
                relativePBC(endToEndDistance[i]);
            }
            double endToEndDistanceModulus = std::sqrt( std::pow(endToEndDistance[0], 2)
                                                      + std::pow(endToEndDistance[1], 2)
                                                      + std::pow(endToEndDistance[2], 2) );
            pOverL += endToEndDistanceModulus/cluster.second.size();
            ++pOverLsamples;
        }
    }
    average_pOverL += pOverL/pOverLsamples;
    outputFile << "Average pOverL at " << simulationTime*simulationTimeStep << ": " << average_pOverL*simulationBoxSide*printingInterval/(simulationTime*simulationTimeStep) << "\n";
}

void IPCsimulation::overrideIPCtypeWithClusterID(std::map<int, std::unordered_set<int>> const& clusters) {
    // overrides the IPC type with a random cluster ID so that it gets printed in the startingstate and trajectory
    int counter = 0;
    for(auto cluster: clusters) {
        ++counter;
        counter %= 24;
        char clusterID = 'A' + counter;
        if (clusterID == 'P')
            clusterID = 'Y';
        if (clusterID == 'Q')
            clusterID = 'Z';
        for (auto ipc: cluster.second) {
            particles[ipc].type = clusterID;
        }
    }
}

void IPCsimulation::printClusterSizes() {
    double norm = printingInterval/simulationTotalDuration;
    std::ofstream clusterSizesFile("siml/averageClusterSizes.out");
    clusterSizesFile << std::scientific << std::setprecision(6);

    int i = 0;
    int integral = 0;
    for(auto n: clusterSizesHistogram) {
        clusterSizesFile << n.first << "\t" << norm*n.second << "\n";
        integral += n.first*n.second;
        ++i;
    }
    clusterSizesFile << "#" << integral << " = " << nIPCs;
    clusterSizesFile.close();
}

void IPCsimulation::computeNematicOrderParameter(std::vector<std::list<int>> const& listOfNeighbours) {
    /* nematic order parameter, I took the definition from
     * Cuetos, van Roij, Dijkstra, Soft Matter, 2008, 4, 757-767
     * https://doi.org/10.1039/B715764A
     */

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

void IPCsimulation::accumulateOrientationsHistogram() {
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
}

void IPCsimulation::printOrientationsHistogram() {
    std::ofstream outputFile("siml/orientationsHistogram.out");
    outputFile << std::scientific << std::setprecision(6);

    double norm = printingInterval/simulationTotalDuration;
    const double binSize = M_PI/orientationHistogramSize;

    for (int p = 0; p < orientationHistogramSize; ++p) {
        for (int a = 0; a < 2*orientationHistogramSize; ++a) {
            outputFile << p*binSize << "\t" << a*binSize << "\t" << orientationsHistogram[p][a]*norm << "\n";
        }
    }
    outputFile.close();
}
