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
    energyTrajectoryFile << "#t\t\t\tT\t\t\tK\t\t\tU\t\t\tE\t\t\trmin\n";

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
    BBinteractionRange = 2*ipcRadius + deltaOverSigma;
    patchRadius = ipcRadius - patchEccentricity;
    inverseTemperature = 1./temperature;

    // output the data for future checks
    printInputFileToOutputFile();

    // scale the lenghts to be in a [0.0:1.0] simulation box
    ipcRadius /= simulationBoxSide;
    ipcDiameter /= simulationBoxSide;
    BBinteractionRange /= simulationBoxSide;
    patchRadius /= simulationBoxSide;
    deltaPotential /= simulationBoxSide;

    // finish processing data
    BBsquaredInteractionRange = std::pow(BBinteractionRange,2);
    BsSquaredInteractionRange = std::pow(ipcRadius + 0.5*deltaPotential + patchRadius,2);
    ssSquaredInteractionRange = std::pow(2.*patchRadius, 2);

    coeff_BB = e_BB / (e_min * deltaPotential);
    coeff_Bs = e_Bs / (e_min * deltaPotential );
    coeff_ss = e_ss / (e_min * deltaPotential);

    // if not restoring, we need to initialize the system here, so that the eccentricities have already been scaled
    if(!stage.inputRestoringPreviousSimulation) {
        outputFile << "Placing " << nIPEs <<  " IPCs on a FCC lattice.\n\n";
        initializeNewConfiguration();
    }

    // cell list compilation
    cells.initialize(1.0, BBinteractionRange, nIPEs);
    outputFile << "Total number of cells: " << cells.getNumberofCells() << std::endl;

    // first computation of the potential
    cells.compileLists(particles);

    for (IPE &ipe: particles) {
        double potential;
        if(computePotentialOfAnIPC(ipe, potential)) {
            std::cerr << "Detected overlap in the initial configuration!\n";
            exit(1);
        }
        ipe.potential = potential;
    }

    computeTotalPotential();
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
    inputFile >> deltaTrans >> deltaRot;
    inputFile.close();
}

//************************************************************************//
void IPEsimulation::printInputFileToOutputFile() {
    outputFile << nIPEs << "\t" << density << "\t" << temperature << "\n";
    outputFile << "\t" << printingInterval << "\t" << simulationTotalDuration << "\n";
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

    for (IPE &ipe: particles) {
        ipe.number = counter++;
        startingConfigurationFile >> unusedPatchName
           >> ipe.cmPosition[0] >> ipe.cmPosition[1] >> ipe.cmPosition[2];
        startingConfigurationFile >> unusedPatchName
           >> ipe.orientation[0] >> ipe.orientation[1] >> ipe.orientation[2];
        startingConfigurationFile >> unusedPatchName
           >> unusedDouble >> unusedDouble >> unusedDouble;
        for (int i: {0, 1, 2}) {
            ipe.orientation[i] = ipe.cmPosition[i] - ipe.orientation[i];
            relativePBC(ipe.orientation[i]);
            ipe.orientation[i] /= patchEccentricity;
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
    for(IPE &ipe: particles) {
        // attempt rotation or translation move
        IPE potentialIPCmove = ipe;
        makeRotationOrTranslationMove(potentialIPCmove, ranGen);

        // check pot diff
        double dU;
        if(computePotentialDifference(potentialIPCmove, dU))
            return;

        // no overlap was detected; if dU negative always accept, otherwise accept with conditional probability
        if (dU <= 0. || ranGen.getRandom01() < std::exp(-dU*inverseTemperature) ) {
            ipe = potentialIPCmove;
            ipe.potential += dU;
        }
    }
}


//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//

void IPEsimulation::makeRotationOrTranslationMove(IPE & ipe, RandomNumberGenerator & ranGen) {
    // choose between a translation and a rotation
    if (ranGen.getRandom55() > 0) {
        // translation
        double delta_x[3];
        generateRandomOrientation(delta_x, ranGen);
        for (int i: {0, 1, 2})
            ipe.cmPosition[i] += deltaTrans*delta_x[i];
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
bool IPEsimulation::computePotentialOfAnIPC(IPE const& ipe, double &dU) {
    dU = 0.;
    int cell = cells.cellNumberFromPosition(ipe);
    const std::list<int> & ipesInCell = cells.getIPCsInCell(cell);
    if (computeInteractionsWithIPEsInTheSameCell(ipe, ipesInCell, dU))
        return true;
    const std::list<int> ipesInNeighbouringCells = cells.getIPCsInNeighbouringCells(cell);
    if (computeInteractionsWithIPEsInNeighbouringCells(ipe, ipesInNeighbouringCells, dU))
        return true;
    return false;
}

//************************************************************************//
bool IPEsimulation::computePotentialDifference(IPE const& ipe, double &dU) {
    // compute interactions of the IPE that was just moved
    double tempativeU = 0.;
    if (computePotentialOfAnIPC(ipe, tempativeU))
        return true;
    dU = tempativeU - ipe.potential;
    return false;
}

//************************************************************************//
bool IPEsimulation::computeInteractionsWithIPEsInTheSameCell(IPE const& ipe, std::list<int> const& ipesInCurrentCell, double& dU) {
    for (auto ins = ipesInCurrentCell.cbegin(); ins != ipesInCurrentCell.cend(); ++ins) {
        if (ipe.number != *ins) { // avoid interaction with the original or itself
            if (computeInteractionsBetweenTwoIPEs(ipe, particles[*ins], dU))
                return true;
        }
    }
    return false;
}

//************************************************************************//
bool IPEsimulation::computeInteractionsWithIPEsInNeighbouringCells(IPE const& ipe, std::list<int> const& ipesInNeighbouringCells, double& dU) {
    for( auto ext = ipesInNeighbouringCells.cbegin(); ext != ipesInNeighbouringCells.cend(); ++ext) {
        if (ipe.number != *ext) { // avoid interaction with the original
                                  // shouldn't be possible, since we don't recompute the list after the move, but just to be safe
            if (computeInteractionsBetweenTwoIPEs(ipe, particles[*ext], dU))
                return true;
        }
    }
    return false;
}

//************************************************************************//
bool IPEsimulation::computeInteractionsBetweenTwoIPEs(const IPE &firstIPE, const IPE &secndIPE, double& dU) {
    double centerCenterSeparation[3];
    for (int i: {0, 1, 2}) {
        centerCenterSeparation[i] = firstIPE.cmPosition[i] - secndIPE.cmPosition[i];
        relativePBC(centerCenterSeparation[i]);
    }
    double centerCenterSeparationModulus = std::pow(centerCenterSeparation[0], 2)
                                         + std::pow(centerCenterSeparation[1], 2)
                                         + std::pow(centerCenterSeparation[2], 2);

    if (centerCenterSeparationModulus < minimumSquaredDistance)
        minimumSquaredDistance = centerCenterSeparationModulus;

    // if the CENTERS are too far, no interactions, skip this couple of IPCs
    if (centerCenterSeparationModulus >= BBsquaredInteractionRange)
        return false;

    // we are inside the interaction range, check the overlap
    if (detectOverlap(firstIPE, secndIPE, centerCenterSeparationModulus))
        return true;

    // no overlap, let's do the real potential computation
    dU += computePotentialBetweenTwoIPEsInsideRange(firstIPE, secndIPE, centerCenterSeparationModulus);

    return false;
}

//************************************************************************//
bool IPEsimulation::detectOverlap(const IPE &firstIPE, const IPE &secndIPE, const double r) {
    if(r < ipcDiameter)
        return true;
    return false;
}

//************************************************************************//
double IPEsimulation::computePotentialBetweenTwoIPEsInsideRange(const IPE &firstIPE, const IPE &secndIPE, const double r) {
    double dU = 0.;
    // compute the interaction between centers
    dU -= coeff_BB*(std::sqrt(r) - BBinteractionRange);

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

            siteSiteSeparationModulus = std::sqrt(siteSiteSeparationModulus);
            dU -= coeff_Bs*(siteSiteSeparationModulus - BBinteractionRange);
        }
        else { // ss
            // if we are too far, no interaction, skip to the next site-site pair
            if (siteSiteSeparationModulus >= ssSquaredInteractionRange)
                continue;

            siteSiteSeparationModulus = std::sqrt(siteSiteSeparationModulus);
            dU -= coeff_ss*(siteSiteSeparationModulus - ssInteractionRange);
        }
    }

    return dU;
}

//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//

void IPEsimulation::computeTotalPotential() {


}

//************************************************************************//
void IPEsimulation::outputSystemTrajectory(std::ofstream & outputTrajectoryFile) {}

//************************************************************************//
void IPEsimulation::outputSystemEnergies(std::ofstream &energyTrajectoryFile) {}
