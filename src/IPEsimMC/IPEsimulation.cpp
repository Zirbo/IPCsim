#include <cstdlib>
#include <iomanip>
#include <iostream>
#include "IPEsimulation.hpp"


IPEsimulation::IPEsimulation(SimulationStage const& stage) {
    // clean up old data and recreate output directory
    if(system("rm -rf siml") != 0) {
        std::cerr << "Unable to delete the old 'siml/' directory with rm -rf. "
                  << "Most likely you have it open somewhere or some program is running in it.\n";
        exit(1);
    }
    if(system("mkdir siml") != 0) {
        std::cerr << "Unable to create a new 'siml/' directory. You'll never see this error message.\n";
        exit(1);
    }

    // open output files
    outputFile.open("siml/output.out");
    energyTrajectoryFile.open("siml/evolution.out");
    energyTrajectoryFile << std::scientific << std::setprecision(6);

    // initialize system
    initializeSystem(stage);

    // print starting configuration and initialize output file
    outputFile << "\nPlot evolution.out to check the evolution of the system.\n";
    outputSystemEnergies(energyTrajectoryFile);


    if (printTrajectoryAndCorrelations) {
        // initialize g(r)
        pairCorrelation.initialize(20, simulationBoxSide, nIPEs);

        // initialize trajectory output file
        trajectoryFile.open("siml/trajectory.xyz");
        trajectoryFile << std::scientific << std::setprecision(24);
        outputSystemTrajectory(trajectoryFile);
    }
}


void IPEsimulation::run() {
    time_t simulationStartTime, simulationEndTime;

    // simulation begins
    time(&simulationStartTime);
    while(simulationTime < simulationTotalDuration) {
        ++simulationTime;
        computeSimulationStep();

        if( simulationTime%printingInterval == 0) {
            // compute and output energies
            outputSystemEnergies(energyTrajectoryFile);
            // output trajectory and compute g(r)
            if (printTrajectoryAndCorrelations) {
                // g(r)
                pairCorrelation.compute(particles);
                outputSystemTrajectory(trajectoryFile);
            }
        }
    }
    // simulation ends
    time(&simulationEndTime);
    outputFile << "The simulation lasted " << difftime (simulationEndTime, simulationStartTime) << " seconds.\n";

    // output final state
    std::ofstream finalStateFile("startingstate.xyz");
    finalStateFile << std::scientific << std::setprecision(24);
    outputSystemTrajectory(finalStateFile);
    finalStateFile.close();

    // output g(r);
    if (printTrajectoryAndCorrelations) {
        const double g_r_integral = pairCorrelation.print("siml/g_r");
        outputFile << "The integral of g(r) is " << g_r_integral << " and is should be equal to the number of particles minus one, " << nIPEs-1 << std::endl;
    }
}


//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//

// Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
void IPEsimulation::generateRandomOrientation(double (&a)[3], RandomNumberGenerator & r) {
    double x,y,quad=2.;
    while ( quad > 1. ) {
        x = r.getRandom11();
        y = r.getRandom11();
        quad = x*x + y*y;
    }
    double norm = 2.*sqrt(1.-quad);  a[0]=x*norm;  a[1]=y*norm;  a[2]=1.-2.*quad;
}

//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//

void IPEsimulation::initializeSystem(SimulationStage const& stage) {

    simulationTime = 0;

    temperature = stage.inputTemperature;
    simulationTotalDuration = stage.inputStageTotalDuration;
    printTrajectoryAndCorrelations = stage.inputPrintTrajectoryAndCorrelations;
    deltaTrans = stage.deltaTrans;
    deltaRot = stage.deltaRot;
    readInputFile();

    // if restoring, read state, so we get access to the real number of IPCs
    if(stage.inputRestoringPreviousSimulation) {
        outputFile << "Resuming a previous simulation. ";
        restorePreviousConfiguration();
        outputFile << "Read " << nIPEs <<  " particles positions and velocities from file.\n";
        // we read nIPCs and simulationBoxSide from the starting configuration, so we can compute the density from them
        density = double(nIPEs)/std::pow(simulationBoxSide, 3);
    }
    else {
        outputFile << "Starting a new simulation.\n";
        // we read nIPCs and density from the input file, so we need to compute the simulationBoxSide from them
        simulationBoxSide = std::cbrt(nIPEs/density);
    }

    // process data
    ipcRadius = 0.5;
    deltaPotential = deltaOverSigma*ipcRadius;
    ipcDiameter = 2.*ipcRadius;
    patchRadius = ipcRadius - patchEccentricity;
    inverseTemperature = 1./temperature;
    BBinteractionRange = 2*ipcRadius + deltaPotential;
    BsinteractionRange = ipcRadius + 0.5*deltaPotential + patchRadius;
    ssInteractionRange = 2.*patchRadius;

    // output the data for future checks
    printInputFileToOutputFile();

    // scale the lenghts to be in a [0.0:1.0] simulation box
    ipcRadius /= simulationBoxSide;
    ipcDiameter /= simulationBoxSide;
    BBinteractionRange /= simulationBoxSide;
    BsinteractionRange /= simulationBoxSide;
    ssInteractionRange /= simulationBoxSide;
    patchEccentricity /= simulationBoxSide;
    patchRadius /= simulationBoxSide;
    deltaPotential /= simulationBoxSide;
    deltaTrans /= simulationBoxSide;

    // finish processing data
    ipcDiameterSquared        = std::pow(ipcDiameter, 2);
    BBsquaredInteractionRange = std::pow(BBinteractionRange,2);
    BsSquaredInteractionRange = std::pow(BsinteractionRange, 2);
    ssSquaredInteractionRange = std::pow(ssInteractionRange, 2);

    coeff_BB = e_BB / (e_min * computeOmega(ipcRadius, ipcRadius, ipcDiameter));
    coeff_Bs = e_Bs / (e_min * computeOmega(patchRadius, ipcRadius, ipcDiameter));
    coeff_ss = e_ss / (e_min * computeOmega(patchRadius, patchRadius, ipcDiameter));


    // if not restoring, we need to initialize the system here, so that the eccentricities have already been scaled
    if(!stage.inputRestoringPreviousSimulation) {
        outputFile << "Placing " << nIPEs <<  " IPCs on a FCC lattice.\n\n";
        initializeNewConfiguration();
    }

    // cell list compilation
    cells.initialize(1.0, BBinteractionRange, nIPEs, cell_lists::SimulationType::MonteCarlo);
    outputFile << "Total number of cells: " << cells.getNumberofCells() << std::endl;

    // first computation of the potential
    cells.compileLists(particles);

    for (IPE &ipe: particles) {
        double potential;
        if(computeFullPotentialOfAnIPE(ipe, potential)) {
            std::cerr << "Detected overlap in the initial configuration!\n";
            exit(1);
        }
        ipe.potential = potential;
    }
}

//************************************************************************//
void IPEsimulation::readInputFile() {
    std::ifstream inputFile("input.in");
    if(inputFile.fail()) {
        std::cerr << "File input.in could not be opened. Aborting.\n";
        exit(1);
    }
    int N1;
    inputFile >> N1 >> density;
    nIPEs = 4*N1*N1*N1;
    inputFile >> printingInterval;
    inputFile >> e_BB >> e_Bs >> e_ss;
    inputFile >> e_min;
    inputFile >> deltaOverSigma >> patchEccentricity;
    inputFile.close();
}

//************************************************************************//
void IPEsimulation::printInputFileToOutputFile() {
    outputFile << nIPEs << "\t" << density << "\t" << temperature << "\n";
    outputFile << printingInterval << "\t" << simulationTotalDuration << "\n";
    outputFile << e_BB << "\t" << e_Bs << "\t" << e_ss << "\n";
    outputFile << e_min << "\n";
    outputFile << deltaOverSigma << "\t" << patchEccentricity << "\n";
    outputFile << deltaTrans << "\t" << deltaRot << "\n";

    outputFile << "\n*****************MC simulation in NVT ensemble for CGDH potential.********************\n";
    outputFile << "\nDensity = " << nIPEs << "/" << std::pow(simulationBoxSide,3) << " = ";
    outputFile << nIPEs/std::pow(simulationBoxSide,3) << " = " << density;
    outputFile << "\nSide = " << simulationBoxSide << ", IPC size in reduced units: " << 1./simulationBoxSide << std::endl;
}

//************************************************************************//
void IPEsimulation::restorePreviousConfiguration() {
    char unusedPatchName;
    double unusedDouble;
    std::ifstream startingConfigurationFile("startingstate.xyz");
    if(startingConfigurationFile.fail()) {
        std::cerr << "File startingstate.xyz could not be opened. Aborting.\n";
        exit(1);
    }
    startingConfigurationFile >> nIPEs >> simulationBoxSide >> unusedDouble;
    nIPEs /= 3;

    particles.resize(nIPEs);
    int counter = 0;

    const double scaledEcc = simulationBoxSide/patchEccentricity;

    for (IPE &ipe: particles) {
        ipe.number = counter++;
        startingConfigurationFile >> unusedPatchName
           >> ipe.cmPosition[0] >> ipe.cmPosition[1] >> ipe.cmPosition[2];
        startingConfigurationFile >> unusedPatchName
           >> ipe.orientation[0] >> ipe.orientation[1] >> ipe.orientation[2];
        startingConfigurationFile >> unusedPatchName
           >> unusedDouble >> unusedDouble >> unusedDouble;
        for (int i: {0, 1, 2}) {
            ipe.orientation[i] -= ipe.cmPosition[i];
            ipe.orientation[i] *= scaledEcc;
        }
    }
    if (counter != nIPEs) {
        std::cerr << "Placed " << counter << " IPCs, expected " << nIPEs << ", quitting.\n";
        exit(1);
    }

    startingConfigurationFile.close();
}

//************************************************************************//
void IPEsimulation::initializeNewConfiguration() {
    particles.resize(nIPEs);
    RandomNumberGenerator rand;

    int N1 = std::cbrt(0.25*nIPEs);
    int N2 = N1*N1;
    int N3 = N2*N1;

    // initialize IPC positions
    for(int i=0;i<N3;i++)
    {
      // FCC is obtained as 4 intersecating SC
        particles[i].number = i;
        particles[i].cmPosition[0] = (i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i].cmPosition[0]);
        particles[i].cmPosition[1] = ((i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i].cmPosition[1]);
        particles[i].cmPosition[2] = (i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(particles[i].cmPosition[2]);

        particles[i+N3].number = i+N3;
        particles[i+N3].cmPosition[0] = (.5 + i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3].cmPosition[0]);
        particles[i+N3].cmPosition[1] = (.5 + (i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3].cmPosition[1]);
        particles[i+N3].cmPosition[2] = (i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(particles[i+N3].cmPosition[2]);

        particles[i+N3+N3].number = i+N3+N3;
        particles[i+N3+N3].cmPosition[0] = (i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3+N3].cmPosition[0]);
        particles[i+N3+N3].cmPosition[1] = (.5 + (i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3+N3].cmPosition[1]);
        particles[i+N3+N3].cmPosition[2] = (.5 + i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(particles[i+N3+N3].cmPosition[2]);

        particles[i+N3+N3+N3].number = i+N3+N3+N3;
        particles[i+N3+N3+N3].cmPosition[0] = (.5 + i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3+N3+N3].cmPosition[0]);
        particles[i+N3+N3+N3].cmPosition[1] = ((i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3+N3+N3].cmPosition[1]);
        particles[i+N3+N3+N3].cmPosition[2] = (.5 + i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(particles[i+N3+N3+N3].cmPosition[2]);

    }
    // initialize patches positions
    for(IPE &ipe: particles) {
        generateRandomOrientation(ipe.orientation, rand);
    }
}

//************************************************************************//
void IPEsimulation::computeSimulationStep() {
    RandomNumberGenerator ranGen;
    minimumSquaredDistance = 1.0;
    for(IPE &ipe: particles) {
        // attempt rotation or translation move
        IPE potentialIPCmove = ipe;
        makeRotationOrTranslationMove(potentialIPCmove, ranGen);

        // check pot diff
        double dU;
        potentialChangeList changes;
        if(computePotentialOfAnIPEmove(potentialIPCmove, dU, changes))
            return;

        // no overlap was detected; if dU negative always accept, otherwise accept with conditional probability
        if (dU <= 0. || ranGen.getRandom01() < std::exp(-dU*inverseTemperature) ) {
            ipe = potentialIPCmove;
            for (potentialChange change: changes) {
                particles[change.particleIndex].potential += change.dU;
            }
        }
    }
}


//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//


bool IPEsimulation::computeFullPotentialOfAnIPE(IPE const& ipe, double &U) {
    U = 0.;
    std::list<int> allNearbyIPEs = findAllTheIPEsInRange(ipe);

    for (auto otherIPE = allNearbyIPEs.cbegin(); otherIPE != allNearbyIPEs.cend(); ++otherIPE) {
        if (ipe.number != *otherIPE) { // avoid interaction with itself
            if (computeInteractionsBetweenTwoIPEs(ipe, particles[*otherIPE], U))
                return true;
        }
    }
    return false;
}

//************************************************************************//
const std::list<int> IPEsimulation::findAllTheIPEsInRange(IPE const& ipe) {
    int cell = cells.cellNumberFromPosition(ipe);
    const std::list<int> ipesInCell(cells.getIPCsInCell(cell));
    const std::list<int> ipesInNeighbouringCells(cells.getIPCsInNeighbouringCells(cell));
    std::list<int> allNearbyIPEs(ipesInCell);
    allNearbyIPEs.insert(allNearbyIPEs.cend(), ipesInNeighbouringCells.cbegin(), ipesInNeighbouringCells.cend());
    return allNearbyIPEs;
}

//************************************************************************//
void IPEsimulation::makeRotationOrTranslationMove(IPE & ipe, RandomNumberGenerator & ranGen) {
    // choose between a translation and a rotation
    if (ranGen.getRandom55() > 0) {
        // translation
        double delta_x[3];
        generateRandomOrientation(delta_x, ranGen);
        for (int i: {0, 1, 2}) {
            ipe.cmPosition[i] += deltaTrans*ranGen.getRandom01()*delta_x[i];
            absolutePBC(ipe.cmPosition[i]);
        }
    } else {
        // rotation
        double delta_n[3];
        generateRandomOrientation(delta_n, ranGen);
        double norm = 0.;
        for (int i: {0, 1, 2}) {
            ipe.orientation[i] += deltaRot*delta_n[i];
            norm += std::pow(ipe.orientation[i], 2);
        }
        norm = std::sqrt(norm);
        for (int i: {0, 1, 2}) {
            ipe.orientation[i] /= norm;
        }
    }
}

//************************************************************************//
bool IPEsimulation::computePotentialOfAnIPEmove(IPE const& ipe, double &dU, potentialChangeList &changes) {
    // compute interactions of the IPE that was just moved
    std::list<int> allNearbyIPEs = findAllTheIPEsInRange(ipe);

    for (auto otherIPE = allNearbyIPEs.cbegin(); otherIPE != allNearbyIPEs.cend(); ++otherIPE) {
        if (ipe.number != *otherIPE) { // avoid interaction with the original or itself
            double Unew = 0.;
            if (computeInteractionsBetweenTwoIPEs(ipe, particles[*otherIPE], Unew))
                return true;

            double Uold = 0.;
            if (computeInteractionsBetweenTwoIPEs(particles[ipe.number], particles[*otherIPE], Uold))
                return true;

            const double deltaUpair = Unew - Uold;
            changes.push_back(potentialChange(*otherIPE, deltaUpair));
            dU += deltaUpair;
        }
    }

    return false;
}

//************************************************************************//
bool IPEsimulation::computeInteractionsWithIPEsInList(const IPE &ipe, std::list<int> const& listOfIPEs, double& dU) {
    for (auto otherIPE = listOfIPEs.cbegin(); otherIPE != listOfIPEs.cend(); ++otherIPE) {
        if (ipe.number != *otherIPE) { // avoid interaction with the original or itself
            if (computeInteractionsBetweenTwoIPEs(ipe, particles[*otherIPE], dU))
                return true;
        }
    }
    return false;
}

//************************************************************************//
bool IPEsimulation::computeInteractionsBetweenTwoIPEs(const IPE &firstIPE, const IPE &secndIPE, double& U) {
    double centerCenterSeparation[3];
    for (int i: {0, 1, 2}) {
        centerCenterSeparation[i] = firstIPE.cmPosition[i] - secndIPE.cmPosition[i];
        relativePBC(centerCenterSeparation[i]);
    }
    double centerCenterSeparationModulus = std::pow(centerCenterSeparation[0], 2)
                                         + std::pow(centerCenterSeparation[1], 2)
                                         + std::pow(centerCenterSeparation[2], 2);

    // if the CENTERS are too far, no interactions, skip this couple of IPCs
    if (centerCenterSeparationModulus >= BBsquaredInteractionRange)
        return false;

    // we are inside the interaction range, check the overlap
    if (detectOverlap(firstIPE, secndIPE, centerCenterSeparationModulus))
        return true;

    // no overlap, let's do the real potential computation
    U += computePotentialBetweenTwoIPEsInsideRange(firstIPE, secndIPE, std::sqrt(centerCenterSeparationModulus));

    if (centerCenterSeparationModulus < minimumSquaredDistance)
        minimumSquaredDistance = centerCenterSeparationModulus;

    return false;
}

//************************************************************************//
bool IPEsimulation::detectOverlap(const IPE &firstIPE, const IPE &secndIPE, const double rSquared) {
    if(rSquared < ipcDiameterSquared)
        return true;
    return false;
}

//************************************************************************//
double IPEsimulation::computePotentialBetweenTwoIPEsInsideRange(const IPE &firstIPE, const IPE &secndIPE, const double r) {
    double U = 0.;
    // compute the interaction between centers
    U += coeff_BB*computeOmega(ipcRadius, ipcRadius, r);

    // compute all the other 8 site-site separations
    double siteSiteSeparation[8][3];
    for (int i: {0, 1, 2}) {
        double fCM = firstIPE.cmPosition[i];
        double fP1 = fCM + patchEccentricity*firstIPE.orientation[i];    absolutePBC(fP1);
        double fP2 = fCM - patchEccentricity*firstIPE.orientation[i];    absolutePBC(fP1);
        double sCM = secndIPE.cmPosition[i];
        double sP1 = sCM + patchEccentricity*secndIPE.orientation[i];    absolutePBC(sP1);
        double sP2 = sCM - patchEccentricity*secndIPE.orientation[i];    absolutePBC(sP1);

        siteSiteSeparation[0][i] = fCM - sP1;
        siteSiteSeparation[1][i] = fCM - sP2;
        siteSiteSeparation[2][i] = fP1 - sCM;
        siteSiteSeparation[3][i] = fP2 - sCM;
        siteSiteSeparation[4][i] = fP1 - sP1;
        siteSiteSeparation[5][i] = fP1 - sP2;
        siteSiteSeparation[6][i] = fP2 - sP1;
        siteSiteSeparation[7][i] = fP2 - sP2;
        for (int j = 0; j < 8; ++j)
            relativePBC(siteSiteSeparation[j][i]);
    }

    for (int j = 0; j < 8; ++j) {
        double siteSiteSeparationModulus = std::pow(siteSiteSeparation[j][0], 2)
                                         + std::pow(siteSiteSeparation[j][1], 2)
                                         + std::pow(siteSiteSeparation[j][2], 2);

        if (j < 4) { // Bs
            // if we are too far, no interaction, skip to the next site-site pair
            if (siteSiteSeparationModulus >= BsSquaredInteractionRange)
                continue;

            U += coeff_Bs*computeOmega(ipcRadius, patchRadius, std::sqrt(siteSiteSeparationModulus));
        }
        else { // ss
            // if we are too far, no interaction, skip to the next site-site pair
            if (siteSiteSeparationModulus >= ssSquaredInteractionRange)
                continue;

            U += coeff_ss*computeOmega(patchRadius, patchRadius, std::sqrt(siteSiteSeparationModulus));
        }
    }

    return U;
}

//************************************************************************//
const double IPEsimulation::computeOmega(const double Ra, const double Rb, const double rab) {
    // BKL paper, formula 18
    if ( rab > Ra+Rb )
        return 0.;
    else if ( rab <= std::fabs(Ra-Rb) )
        return 8.*std::pow(std::min(Ra,Rb),3);
    else {
        const double tempSum = (Ra*Ra-Rb*Rb)/(2.*rab);
        return 2.*( (2.*Ra+tempSum+rab/2.)*std::pow(Ra-tempSum-rab/2.,2)
                  + (2.*Rb-tempSum+rab/2.)*std::pow(Rb+tempSum-rab/2.,2) );
    }
}

//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//

void IPEsimulation::outputSystemTrajectory(std::ofstream & outputTrajectoryFile) {
    outputTrajectoryFile << 3*nIPEs << "\n" << simulationBoxSide << "\t" << simulationTime;
    for (IPE &ipe: particles) {
        outputTrajectoryFile << "\n" << "C" << "\t"
                             << ipe.cmPosition[0] << "\t" << ipe.cmPosition[1] << "\t" << ipe.cmPosition[2];
        outputTrajectoryFile << "\n" << "P" << "\t"
                             << ipe.cmPosition[0] + patchEccentricity*ipe.orientation[0] << "\t"
                             << ipe.cmPosition[1] + patchEccentricity*ipe.orientation[1] << "\t"
                             << ipe.cmPosition[2] + patchEccentricity*ipe.orientation[2];
        outputTrajectoryFile << "\n" << "Q" << "\t"
                             << ipe.cmPosition[0] - patchEccentricity*ipe.orientation[0] << "\t"
                             << ipe.cmPosition[1] - patchEccentricity*ipe.orientation[1] << "\t"
                             << ipe.cmPosition[2] - patchEccentricity*ipe.orientation[2];
    }
    outputTrajectoryFile << std::endl;
}

//************************************************************************//
void IPEsimulation::outputSystemEnergies(std::ofstream &energyTrajectoryFile) {
    double U = 0.;
    for (IPE const& ipe: particles) {
        U += ipe.potential;
    }
    U /= 2.*nIPEs;

    energyTrajectoryFile << simulationTime << "\t" << U <<"\t" << simulationBoxSide*std::sqrt(minimumSquaredDistance) << std::endl;
}
