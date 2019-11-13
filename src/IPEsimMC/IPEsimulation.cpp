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
        pairCorrelation.initialize(20, simulationBoxSide, nIPCs);

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
        outputFile << "The integral of g(r) is " << g_r_integral << " and is should be equal to the number of particles minus one, " << nIPCs-1 << std::endl;
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

// compute overlap values between two spheres of radius Ra and Rb, at distance rab
double IPEsimulation::computeOmega(double Ra, double Rb, double rab) {
    // BKL paper, formula 18
    if ( rab > Ra+Rb )
        return 0.;
    else if ( rab <= std::fabs(Ra-Rb) )
        return 8.*std::pow(std::min(Ra,Rb),3);
    else {
        const double tempSum = (Ra*Ra-Rb*Rb)/(2.*rab);
        return 2.*( (2.*Ra+tempSum+rab/2.)*pow(Ra-tempSum-rab/2.,2)
                  + (2.*Rb-tempSum+rab/2.)*pow(Rb+tempSum-rab/2.,2) );
    }
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
        outputFile << "Read " << nIPCs <<  " particles positions and velocities from file.\n";
        // we read nIPCs and simulationBoxSide from the starting configuration, so we can compute the density from them
        density = double(nIPCs)/std::pow(simulationBoxSide, 3);
    }
    else {
        outputFile << "Starting a new simulation.\n";
        // we read nIPCs and density from the input file, so we need to compute the simulationBoxSide from them
        simulationBoxSide = std::cbrt(nIPCs/density);
    }

    // process data
    ipcRadius = firstPatchEccentricity + firstPatchRadius;  // works for both 2patch and Janus
    interactionRange = 2*ipcRadius;


    // output the data for future checks
    outputFile << nIPCs << "\t" << density << "\t" << temperature << "\n";
    outputFile << "\t" << printingInterval << "\t" << simulationTotalDuration << "\n";
    outputFile << e_BB << "\t" << e_Bs1 << "\t" << e_Bs2 << "\n";
    outputFile << e_s1s2 << "\t" << e_s1s1 << "\t" << e_s2s2 << "\n";
    outputFile << e_min << "\n";
    outputFile << firstPatchEccentricity << "\t" << firstPatchRadius << "\n";
    outputFile << secndPatchEccentricity << "\t" << secndPatchRadius << "\n";
    outputFile << forceAndEnergySamplingStep << "\n";

    outputFile << "\n*****************MC simulation in NVT ensemble for CGDH potential.********************\n";
    outputFile << "\nDensity = " << nIPCs << "/" << std::pow(simulationBoxSide,3) << " = ";
    outputFile << nIPCs/std::pow(simulationBoxSide,3) << " = " << density;
    outputFile << "\nSide = " << simulationBoxSide << ", IPC size in reduced units: " << 1./simulationBoxSide << std::endl;

    // scale the lenghts to be in a [0.0:1.0] simulation box
    ipcRadius /= simulationBoxSide;
    interactionRange /= simulationBoxSide;
    firstPatchRadius /= simulationBoxSide;
    firstPatchEccentricity /= simulationBoxSide;
    secndPatchRadius /= simulationBoxSide;
    secndPatchEccentricity /= simulationBoxSide;
    forceAndEnergySamplingStep /= simulationBoxSide;

    // finish processing data
    squaredInteractionRange = std::pow(interactionRange,2);

    // if not restoring, we need to initialize the system here, so that the eccentricities have already been scaled
    if(!stage.inputRestoringPreviousSimulation) {
        outputFile << "Placing " << nIPCs <<  " IPCs on a FCC lattice.\n\n";
        initializeNewConfiguration();
    }

    // cell list compilation
    cells.initialize(1., interactionRange, nIPCs);
    outputFile << "Total number of cells: " << cells.getNumberofCells() << std::endl;

    // first computation of the potential
    cells.compileLists(particles);
    computePotential();
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
    nIPCs = 4*N1*N1*N1;
    inputFile >> printingInterval;
    inputFile >> e_BB >> e_Bs1 >> e_Bs2;
    inputFile >> e_s1s1 >> e_s2s2 >> e_s1s2;
    inputFile >> e_min;
    inputFile >> firstPatchEccentricity >> firstPatchRadius;
    inputFile >> secndPatchEccentricity >> secndPatchRadius;
    inputFile >> forceAndEnergySamplingStep;
    inputFile.close();

    // patch geometry integrity check
    if ( std::abs( (firstPatchEccentricity+firstPatchRadius)-(secndPatchEccentricity+secndPatchRadius) ) >= 1e-10 ) {
        std::cerr << firstPatchEccentricity << "+" << firstPatchRadius << "=" << firstPatchEccentricity+firstPatchRadius << "-";
        std::cerr << secndPatchEccentricity << "+" << secndPatchRadius << "=" << secndPatchEccentricity+secndPatchRadius << "=\n";
        std::cerr << (firstPatchEccentricity+firstPatchRadius)-(secndPatchEccentricity+secndPatchRadius) << std::endl;
        std::cerr << "eccentricities and radii are not consistent!\n";
        exit(1);
    }
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
    startingConfigurationFile >> nIPCs >> simulationBoxSide >> unusedDouble;
    nIPCs /= 3;

    particles.resize(nIPCs);
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
            ipe.orientation[i] /= firstPatchEccentricity;
        }
    }
    if (counter != nIPCs) {
        std::cerr << "Placed " << counter << " IPCs, expected " << nIPCs << ", quitting.\n";
        exit(1);
    }

    startingConfigurationFile.close();
}

//************************************************************************//
void IPEsimulation::initializeNewConfiguration() {
    particles.resize(nIPCs);
    RandomNumberGenerator rand;

    int N1 = std::cbrt(0.25*nIPCs);
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
        IPE backup = ipe;
        makeRotationOrTranslationMove(ipe, ranGen);

        // check pot diff
        double dU = computePotentialDifference(ipe);
        // if dU < 0 we are done, otherwise revert.
        if (dU > 0.) {
            ipe = backup;
        }
    }
}


//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//

double IPEsimulation::makeRotationOrTranslationMove(IPE & ipe, RandomNumberGenerator & ranGen) {
    // choose between a rotation and a translation
    if (ranGen.getRandom55() > 0) {
        // rotational move

    } else {
        // translation

    }
}
//************************************************************************//
double IPEsimulation::computePotentialDifference(IPE const& ipe) {
  /*  int cell = cells.cellNumberFromPosition(ipe);
    const std::list<int> & ipesInCell = cells.getIPCsInCell(cell);
    const std::list<int> ipesInNeighbouringCells = cells.getIPCsInNeighbouringCells(cell);
    double dU = 0.;
    dU += computeInteractionsWithIPEsInTheSameCell(ipe, ipesInCell);
    dU += computeInteractionsWithIPEsInNeighbouringCells(ipe, ipesInNeighbouringCells);
    return dU;*/
    return 2;
}

//************************************************************************//
double IPEsimulation::computeInteractionsWithIPEsInTheSameCell(IPE const& ipe, std::list<int> const& ipesInCurrentCell) {
    return 2;
}

//************************************************************************//
double IPEsimulation::computeInteractionsWithIPEsInNeighbouringCells(IPE const& ipe, std::list<int> const& ipesInNeighbouringCells) {
    return 2;
}

//************************************************************************//
double IPEsimulation::computeInteractionsBetweenTwoIPEs(const int firstIPC, const int secndIPC) {
    return 2;
}


//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//

void IPEsimulation::computePotential() {}

//************************************************************************//
void IPEsimulation::outputSystemTrajectory(std::ofstream & outputTrajectoryFile) {}

//************************************************************************//
void IPEsimulation::outputSystemEnergies(std::ofstream &energyTrajectoryFile) {}
