#include "IPCpostprocess.hpp"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>


IPCpostprocess::IPCpostprocess(const int inputNumberOfPatches, const int inputNumberOfSubSimulations, std::string const& directoryName) {
    // input checks
    if (inputNumberOfPatches != 1 && inputNumberOfPatches != 2) {
        std::cerr << "At the moment only 1 and 2 patches are supported.\n";
        exit(1);
    }
    numberOfPatches = inputNumberOfPatches;
    if (inputNumberOfSubSimulations <= 0) {
        std::cerr << "The number of subsimulations needs to be a positive integer.\n";
        exit(1);
    }
    numberOfSubSimulations = inputNumberOfSubSimulations;

    // clean up old data and recreate output directory
    if(system("rm -rf analysis") != 0) {
        std::cerr << "Unable to delete the old 'analysis/' directory with rm -rf. "
                  << "Most likely you have it open somewhere or some program is running in it.\n";
        exit(1);
    }
    if(system("mkdir analysis") != 0) {
        std::cerr << "Unable to create a new 'analysis/' directory. You'll never see this error message.\n";
        exit(1);
    }

    // open simulation output file and trajectory
    const std::string trajectoryFileName = directoryName + "/trajectory.xyz";
    trajectoryFile.open(trajectoryFileName);
    if(trajectoryFile.fail()) {
        std::cerr << "File " << trajectoryFileName << " could not be opened. Aborting.\n";
        exit(1);
    }

    readOutputFile(directoryName);

    // initialize the containers
    ipcCentersPreviousPositions.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCentersCurrentPositions.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcInitialOrientations.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCurrentOrientations.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCentersInitialVelocities.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCentersCurrentVelocities.resize(nIPCs, {0.0, 0.0, 0.0});

    displacementOfEachIPCs.resize(nIPCs, {0.0, 0.0, 0.0});
    meanSquaredDisplacement.resize(simulationDurationInIterations, 0.);
    orientationAutocorrelation.resize(subSimulationDuration, 0.);
    velocityAutocorrelation.resize(subSimulationDuration, 0.);
    orientationHistogramSize = 40;
    totalCollectedOrientations = 0;
    orientationHistogram3D.resize(orientationHistogramSize, std::vector<unsigned long>(orientationHistogramSize));
    orientationHistogram2D.resize(orientationHistogramSize, std::vector<unsigned long>(orientationHistogramSize));
}

void IPCpostprocess::run() {
    for (int subSym = 0; subSym < numberOfSubSimulations; ++subSym) {
        for (int subSymSnapshot = 0; subSymSnapshot < subSimulationDuration; ++subSymSnapshot) {
            const int absoluteSnapshot = subSym*subSimulationDuration + subSymSnapshot;
            readSnapshot(absoluteSnapshot);
            if (subSymSnapshot == 0)
                updateInitialOrientationAndVelocites();
            if (absoluteSnapshot == 0) {
                updatePreviousPositions(); // for the initial computation of the MSD...
                computeStaticProperties();
            }
            // compute dynamic properties, or static properties that we average over the entire simulation run
            computeMSD(absoluteSnapshot);
            computeAutocorrelations(subSymSnapshot);
            accumulateTypicalOrientations();
            updatePreviousPositions();
        }
        std::cout << "Subsym " << subSym+1 << " of " << numberOfSubSimulations << " finished\n";
    }

    printAutocorrelations();
    printMSD();
    printTypicalOrientations();

    trajectoryFile.close();
}


void IPCpostprocess::readSnapshot(const int snapshotNumber) {
    runConsistencyChecks(snapshotNumber);

    for(int i = 0; i < nIPCs; ++i) {
        char uselessChar;
        double useless;
        // IPC center
        trajectoryFile >> uselessChar
                >> ipcCentersCurrentPositions[i][0] >> ipcCentersCurrentPositions[i][1] >> ipcCentersCurrentPositions[i][2]
                >> ipcCentersCurrentVelocities[i][0] >> ipcCentersCurrentVelocities[i][1] >> ipcCentersCurrentVelocities[i][2];
        // first patch
        trajectoryFile >> uselessChar
                >> ipcCurrentOrientations[i][0] >> ipcCurrentOrientations[i][1] >> ipcCurrentOrientations[i][2]
                >> useless >> useless >> useless;
        // second patch
        if (numberOfPatches == 2) {
            trajectoryFile >> uselessChar;
            // ignore the rest of the line
            trajectoryFile.ignore(500, '\n');
        }
        for (int j: {0, 1, 2}) {
            // compute orientations
            ipcCurrentOrientations[i][j] -= ipcCentersCurrentPositions[i][j];
            relativePBC(ipcCurrentOrientations[i][j]);
            ipcCurrentOrientations[i][j] /= firstPatchEccentricity;
        }
    }
}

void IPCpostprocess::updateInitialOrientationAndVelocites() {
    // set the current as the initial state
    for(int i = 0; i < nIPCs; ++i) {
        for (int j: {0, 1, 2}) {
            ipcInitialOrientations[i][j] = ipcCurrentOrientations[i][j];
            ipcCentersInitialVelocities[i][j] = ipcCentersCurrentVelocities[i][j];
        }
    }
}

void IPCpostprocess::updatePreviousPositions() {
    // set the current as the previous iteration state
    for(int i = 0; i < nIPCs; ++i) {
        for (int j: {0, 1, 2}) {
            ipcCentersPreviousPositions[i][j] = ipcCentersCurrentPositions[i][j];
        }
    }

}

void IPCpostprocess::runConsistencyChecks(const int snapshotNumber) {
    int nIPCsCheck;
    double simulationBoxSideCheck;
    double time;
    trajectoryFile >> nIPCsCheck >> simulationBoxSideCheck >> time;
    if ((1+numberOfPatches)*nIPCs != nIPCsCheck) {
        std::cerr << "Inconsistency in the number of particles.\n"
                  << (1+numberOfPatches)*nIPCs << " != " << nIPCsCheck << ".\n";
        exit(1);
    }
    if ( std::fabs(simulationBoxSide - simulationBoxSideCheck) > 1e-2) {
        std::cerr << "Inconsistency in the simulation box side.\n"
                  << simulationBoxSide << " != " << simulationBoxSideCheck << ".\n";
        exit(1);
    }
    if (snapshotNumber*printingInterval != time) {
        std::cerr << "Inconsistency in the time.\n"
                  << snapshotNumber*printingInterval << " != " << time << ".\n";
        exit(1);
    }
}

void IPCpostprocess::readOutputFile(std::string const& directoryName) {
    // open file
    const std::string outputFileName = directoryName + "/output.out";
    std::ifstream simulationOutputFile(outputFileName);
    if(simulationOutputFile.fail()) {
        std::cerr << "File " << outputFileName << " could not be opened. Aborting.\n";
        exit(1);
    }

    // read simulation parameters from the simulation output file
    simulationOutputFile.ignore(200, '\n');
    double useless;
    simulationOutputFile >> nIPCs >> density >> temperature
            >> simulationTimeStep >> printingInterval >> simulationTotalDuration
            >> useless >> useless >> useless
            >> useless >> useless >> useless
            >> useless
            >> firstPatchEccentricity >> firstPatchRadius
            >> secndPatchEccentricity >> secndPatchRadius;

    //close file
    simulationOutputFile.close();

    // process data
    simulationBoxSide = std::cbrt(nIPCs/density);
    simulationDurationInIterations = (int)simulationTotalDuration/printingInterval;
    subSimulationDuration = (int)simulationDurationInIterations/numberOfSubSimulations;

    // scale lenghts
    firstPatchEccentricity /= simulationBoxSide;
    firstPatchRadius /= simulationBoxSide;

    // compute working parameters
    squaredInteractionRange = std::pow(2*(firstPatchEccentricity + firstPatchRadius), 2);
}

void IPCpostprocess::computeMSD(const int snapshotNumber) {
    for(int i = 0; i < nIPCs; ++i) {
        for (int j: {0, 1, 2}) {
            double delta_xj = ipcCentersCurrentPositions[i][j] - ipcCentersPreviousPositions[i][j];
            relativePBC(delta_xj);
            displacementOfEachIPCs[i][j] += delta_xj;
            meanSquaredDisplacement[snapshotNumber] += std::pow(displacementOfEachIPCs[i][j],2);
        }
    }
    meanSquaredDisplacement[snapshotNumber] /= nIPCs;
}

void IPCpostprocess::printMSD() {
    std::ofstream meanSquaredDisplFile("analysis/meanSquaredDisplacement.out");
    meanSquaredDisplFile << std::scientific << std::setprecision(6);
    for(int snapshotNumber = 0; snapshotNumber < simulationDurationInIterations; ++snapshotNumber) {
        meanSquaredDisplFile << snapshotNumber*printingInterval << "\t" << meanSquaredDisplacement[snapshotNumber] << "\n";
    }
    meanSquaredDisplFile.close();
}

void IPCpostprocess::computeAutocorrelations(const int snapshotNumber) {
    for(int i = 0; i < nIPCs; ++i) {
        for (int j: {0, 1, 2}) {
            // compute autocorrelations
            orientationAutocorrelation[snapshotNumber] += ipcCurrentOrientations[i][j]*ipcInitialOrientations[i][j];
            velocityAutocorrelation[snapshotNumber] += ipcCentersCurrentVelocities[i][j]*ipcCentersInitialVelocities[i][j];
        }
    }
}

void IPCpostprocess::printAutocorrelations() {
    std::ofstream autocorrelationsFile("analysis/autocorrelations.out");
    autocorrelationsFile << std::scientific << std::setprecision(6);

    const double orientationScalingFactor = 1./nIPCs/numberOfSubSimulations;
    const double velocityScalingFactor = 1./velocityAutocorrelation[0];
    for (int i = 0; i < subSimulationDuration; ++i) {
        autocorrelationsFile << i*printingInterval << "\t"
                             << orientationAutocorrelation[i]*orientationScalingFactor << "\t"
                             << velocityAutocorrelation[i]*velocityScalingFactor << "\n";
    }
    autocorrelationsFile.close();
}

void IPCpostprocess::accumulateTypicalOrientations() {
    const double azimuthConversionFactor = double(orientationHistogramSize)/M_PI;
    const double polarConversionFactor = double(orientationHistogramSize)/(2.*M_PI);
    const int halfHistogramSize = int(orientationHistogramSize/2);

    for (auto const& ipcOrientation: ipcCurrentOrientations) {
        ++totalCollectedOrientations;

        const int xBin = halfHistogramSize + (int)halfHistogramSize*ipcOrientation[0];
        const int yBin = halfHistogramSize + (int)halfHistogramSize*ipcOrientation[1];
        ++orientationHistogram2D[xBin][yBin];

        double azimuthAngle = std::acos(ipcOrientation[2]);
        double polarAngle = std::atan2(ipcOrientation[1], ipcOrientation[0]);
        if (polarAngle < 0.)
            polarAngle += M_PI;
        const int azimuthAngleBin = (int) std::floor(azimuthAngle*azimuthConversionFactor);
        const int polarAngleBin = (int) std::floor(polarAngle*polarConversionFactor);
        ++orientationHistogram3D[polarAngleBin][azimuthAngleBin];
    }
}

void IPCpostprocess::printTypicalOrientations() {
    std::ofstream typicalOrientationsFile3D("analysis/typicalOrientations3D.out");
    typicalOrientationsFile3D << std::scientific << std::setprecision(6);

    std::ofstream typicalOrientationsFile2D("analysis/typicalOrientations2D.out");
    typicalOrientationsFile2D << std::scientific << std::setprecision(6);

    if (simulationDurationInIterations*nIPCs != totalCollectedOrientations) {
        std::cerr << "Consistency check failed.\n" << simulationDurationInIterations << "x" << nIPCs << " != " << totalCollectedOrientations << "!\n";
        exit(1);
    }
    const double azimuthConversionFactor = M_PI/orientationHistogramSize;
    const double polarConversionFactor = 2.*M_PI/orientationHistogramSize;
    const int halfHistogramSize = (int)orientationHistogramSize/2;
    const double inverseHalfHistogramSize = 1./halfHistogramSize;
    const double inverseHistogramSize = 1./orientationHistogramSize;
    const double inverseTotalCollectedOrientations = 1./totalCollectedOrientations;
    for (int azimuthBin = 0; azimuthBin <= orientationHistogramSize; ++azimuthBin) {
        for (int polarBin = 0; polarBin <= orientationHistogramSize; ++polarBin) {
            const double azimuthAngle = azimuthBin*azimuthConversionFactor;
            const double polarAngle = polarBin*polarConversionFactor;
            // correction for the full angles that are needed by gnuplot!
            const int printingAzimuthBin = (azimuthBin == orientationHistogramSize)? 0 : azimuthBin;
            const int printingPolarBin = (polarBin == orientationHistogramSize)? 0 : polarBin;
            typicalOrientationsFile3D << azimuthAngle << "\t" << polarAngle << "\t" << inverseTotalCollectedOrientations*orientationHistogram3D[printingAzimuthBin][printingPolarBin] << "\n";

            if (azimuthBin == orientationHistogramSize || polarBin == orientationHistogramSize)
                continue;
            const double x = (azimuthBin-halfHistogramSize)*inverseHalfHistogramSize + inverseHistogramSize;
            const double y = (polarBin-halfHistogramSize)*inverseHalfHistogramSize + inverseHistogramSize;
            const double z = inverseTotalCollectedOrientations*orientationHistogram2D[azimuthBin][polarBin];
            typicalOrientationsFile2D << x << "\t" << y << "\t" << z << "\t" << std::log(z) << "\n";
        }
        typicalOrientationsFile3D << std::endl;
    }
    typicalOrientationsFile3D.close();
    typicalOrientationsFile2D.close();
}

void IPCpostprocess::computeStaticProperties() {
    const std::vector<std::list<int>> listOfNeighbours = computeListOfNeighbours();
    computeAndPrintHistogramOfNeighbours(listOfNeighbours);
    computeNematicOrderParameter(listOfNeighbours);
    computePcosTheta(listOfNeighbours);
    doClusterAnalysis(listOfNeighbours);
}

std::vector<std::list<int>> IPCpostprocess::computeListOfNeighbours() {
    std::vector<std::list<int>> listOfNeighbours(nIPCs);
    for (int i = 0; i < nIPCs; ++i) {
        for (int j = i + 1; j < nIPCs; ++j) {
            double distance = 0.;
            for (int d: {0, 1, 2}) {
                double distance_d = ipcCentersPreviousPositions[i][d] - ipcCentersPreviousPositions[j][d];
                relativePBC(distance_d);
                distance += std::pow(distance_d, 2);
            }
            // if they are too far, it does not count
            if (distance > squaredInteractionRange)
                continue;
            //else
            listOfNeighbours[i].push_back(j);
            listOfNeighbours[j].push_back(i);
        }
    }
    return listOfNeighbours;
}

void IPCpostprocess::computeAndPrintHistogramOfNeighbours(std::vector<std::list<int>> const& listOfNeighbours) {
    // compute how many neighbours each particle has
    std::vector<int> numberOfNeighbours(nIPCs, 0);
    for (int i = 0; i < nIPCs; ++i)
        numberOfNeighbours[i] = listOfNeighbours[i].size();
    // compute the histogram of neighbours
    const int maxNumberOfNeighbours = *std::max_element(numberOfNeighbours.cbegin(), numberOfNeighbours.cend());
    std::vector<int> histogramOfNeighbours(maxNumberOfNeighbours+1, 0);
    for(int neighboursOfThisParticle: numberOfNeighbours)
        ++histogramOfNeighbours[neighboursOfThisParticle];

    // and print it
    std::ofstream numberOfNeighboursFile("analysis/numberOfNeighbours.out");
    numberOfNeighboursFile << std::scientific << std::setprecision(6);
    for(int i = 0; i < maxNumberOfNeighbours;  ++i) {
        numberOfNeighboursFile << i << "\t" << histogramOfNeighbours[i] << "\n";
    }
    numberOfNeighboursFile.close();
}

void IPCpostprocess::computeNematicOrderParameter(std::vector<std::list<int>> const& listOfNeighbours) {
    std::ofstream nematicOrderParameterFile("analysis/nematicOrderParameter.out");
    nematicOrderParameterFile << std::scientific << std::setprecision(6);

    double globalAverage = 0.;
    // loop on the lists of neighbours
    for (int i = 0; i < nIPCs; ++i) {
        double modulusNOPi = 0.;
        // loop on the neighbours inside the list
        for (int j: listOfNeighbours[i]) {
            double modulusNOPij = 0.;
            for (int d: {0, 1, 2})
                modulusNOPij += ipcCurrentOrientations[i][d]*ipcCurrentOrientations[j][d];
            modulusNOPi += std::pow(modulusNOPij,2);
        }
        // normalize and print
        if(!listOfNeighbours[i].empty())
            modulusNOPi = 1.5*modulusNOPi/listOfNeighbours[i].size() - 0.5;
        nematicOrderParameterFile << i << "\t" << modulusNOPi << "\n";
        globalAverage += modulusNOPi;
    }
    nematicOrderParameterFile.close();
    globalAverage /= nIPCs;
    std::cout << "Global average of the nematic order parameter: " << globalAverage << "!\n";
}

void IPCpostprocess::computePcosTheta(std::vector<std::list<int>> const& listOfNeighbours) {
    std::ofstream nematicOrderParameterFile("analysis/PcosTheta.out");
    nematicOrderParameterFile << std::scientific << std::setprecision(6);

    double globalAverage = 0.;
    double globalAverageAbsoluteValue = 0.;
    // loop on the lists of neighbours
    for (int i = 0; i < nIPCs; ++i) {
        double modulusPcosTheta = 0.;
        double modulusPcosThetaAbsoluteValue = 0.;
        // loop on the neighbours inside the list
        for (int j: listOfNeighbours[i]) {
            double modulusPcosTheta_ij = 0.;
            for (int d: {0, 1, 2})
                modulusPcosTheta_ij += ipcCurrentOrientations[i][d]*ipcCurrentOrientations[j][d];
            modulusPcosTheta += modulusPcosTheta_ij;
            modulusPcosThetaAbsoluteValue += std::abs(modulusPcosTheta_ij);

        }
        // normalize and print
        if(!listOfNeighbours[i].empty()) {
            modulusPcosTheta = modulusPcosTheta/listOfNeighbours[i].size();
            modulusPcosThetaAbsoluteValue = modulusPcosThetaAbsoluteValue/listOfNeighbours[i].size();
        }
        nematicOrderParameterFile << i << "\t" << modulusPcosTheta << "\t" << modulusPcosThetaAbsoluteValue << "\n";
        globalAverage += modulusPcosTheta;
        globalAverageAbsoluteValue += modulusPcosThetaAbsoluteValue;
    }
    nematicOrderParameterFile.close();
    globalAverage /= nIPCs;
    globalAverageAbsoluteValue /= nIPCs;
    std::cout << "Global average of PcosTheta: " << globalAverage << "!\n";
    std::cout << "Global average of |PcosTheta|: " << globalAverageAbsoluteValue << "!\n";
}

void IPCpostprocess::doClusterAnalysis(std::vector<std::list<int>> const& listOfNeighbours) {
    // construct a matrix with the particles parallel to each other...
 /*   std::vector<std::vector<bool>> isParallel(nIPCs, std::vector<bool>(nIPCs, false));
    for (int i = 0; i < nIPCs; ++i) {
        for (int j: listOfNeighbours[i]) {
            double scalarProduct_ij = 0.;
            for (int d: {0, 1, 2})
                scalarProduct_ij += ipcCurrentOrientations[i][d]*ipcCurrentOrientations[j][d];
            scalarProduct_ij = std::pow(scalarProduct_ij,2);
            if(scalarProduct_ij > .8)
                isParallel[i][j] = true;
        }
    }*/
    // e mo che cazzo ci volevo fare?
}
