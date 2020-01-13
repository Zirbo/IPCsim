#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "IPCsimulation.hpp"


//************************************************************************//
IPCsimulation::IPCsimulation(SimulationStage const& stage) {
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
    energyTrajectoryFile << std::scientific << std::setprecision(10);
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
        outputSystemTrajectory(trajectoryFile, printForces);

        // initialize mean square displacement file and helpers
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
}

//************************************************************************//
double IPCsimulation::run() {
    time_t simulationStartTime, simulationEndTime;

    const size_t simulationDurationInIterations = (size_t)simulationTotalDuration/simulationTimeStep;
    const double printingIntervalDouble = printingInterval/simulationTimeStep;
    const size_t printingIntervalInIterations = (size_t)printingIntervalDouble;

    double averageTemperature = 0.;
    double averageSquaredTemperature = 0.;
    double averagePotentialEnergy = 0.;
    int sampleSizeForAverages = 0;

    // simulation begins
    time(&simulationStartTime);
    while(simulationTime < simulationDurationInIterations) {
        ++simulationTime;
        computeTrajectoryStep();

        if( simulationTime%printingIntervalInIterations == 0) {
            // compute and output energies
            computeSystemEnergy();
            outputSystemEnergies(energyTrajectoryFile);

            // output trajectory and compute analysis
            if (printTrajectoryAndCorrelations) {
                pairCorrelation.compute(particles);
                outputSystemTrajectory(trajectoryFile, printForces);
                computeMSD();
                computeStaticProperties();
            }
            // compute averages
            ++sampleSizeForAverages;
            averageTemperature += temperature;
            averageSquaredTemperature += temperature*temperature;
            averagePotentialEnergy += potentialEnergy;
        }
    }
    // simulation ends
    time(&simulationEndTime);
    outputFile << "The simulation lasted " << difftime (simulationEndTime, simulationStartTime) << " seconds.\n";

    // output final state
    std::ofstream finalStateFile("startingstate.xyz");
    finalStateFile << std::scientific << std::setprecision(24);
    outputSystemTrajectory(finalStateFile, false);
    finalStateFile.close();

    // check that total momentum is still zero and print final stuff
    double pcm [3];
    computeSystemMomentum(pcm);
    outputFile << "Residual momentum of the whole system = ( " << pcm[0]*simulationBoxSide << ", " << pcm[1]*simulationBoxSide << ", " << pcm[2]*simulationBoxSide << " ).\n" << std::endl;
    averageTemperature /= sampleSizeForAverages;
    averagePotentialEnergy /= sampleSizeForAverages;
    averageSquaredTemperature /= sampleSizeForAverages;
    double temperatureVariance = std::sqrt(averageSquaredTemperature - std::pow(averageTemperature,2));
    outputFile << "Average kT during the simulation run = " << averageTemperature << std::endl;
    outputFile << "Standard deviation of kT during the simulation run = " << std::sqrt(temperatureVariance) << std::endl;
    outputFile << "Average potential energy during the simulation run = " << averagePotentialEnergy/nIPCs << std::endl;

    // output analysis
    if (printTrajectoryAndCorrelations) {
        printStaticProperties();
    }

    return averageTemperature;
}

//************************************************************************//
void IPCsimulation::printPotentials() {
    // clean up unneeded shit
    outputFile.close();
    trajectoryFile.close();
    energyTrajectoryFile.close();
    if(system("rm -rf siml") != 0) {
        std::cerr << "Unable to delete the just created 'siml/' directory with rm -rf. "
                  << "Most likely you have it open somewhere or some program is running in it.\n";
        exit(1);
    }
    if(system("rm -rf potentials_for_lammps potentials") != 0) {
        std::cerr << "Unable to delete the old 'potentials_for_lammps' and/or 'potentials' directory with rm -rf. "
                  << "Most likely you have it open somewhere, or some program is running in it.\n";
        exit(1);
    }

    int choice;
    std::cout << "Which potentials do you want to print?\n 1 Silvano\n 2 Emanuela\n 3 LAMMPS\n";
    std::cin >> choice;


    if (choice == 1) {
        // RAW SITE-SITE
        if(system("mkdir potentials") != 0) {
            std::cerr << "Unable to create a new 'potentials_raw_sitesite/' directory. You'll never see this error message.\n";
            exit(1);
        }

        int potentialPrintingStep = 10000;
        std::cout << "You chose raw site-site potentials.\n";
        printRawSiteSitePotentials(potentialPrintingStep);;
    }
    else if (choice == 2) {
        // CONTOUR PLOTS
        if(system("mkdir potentials") != 0) {
            std::cerr << "Unable to create a new 'potentials/' directory. You'll never see this error message.\n";
            exit(1);
        }

        int potentialPrintingStep = 10000;
        std::cout << "You chose contour plots\n";
        printPotentialsToFileForVisualization(potentialPrintingStep);;
    }
    else if (choice == 3) {
        // LAMMPS
        if(system("mkdir potentials_for_lammps") != 0) {
            std::cerr << "Unable to create a new 'potentials_for_lammps/' directory. You'll never see this error message.\n";
            exit(1);
        }
        int potentialPrintingStep = 20000;
        int cutoffValue = 100;
        std::cout << "You chose LAMMPS\n"
                  << "Your potential is defined every " << forceAndEnergySamplingStep*simulationBoxSide
                  << " and until " << interactionRange*simulationBoxSide << ".\n"
                  << "How often do you want to print, in integer multiples of "
                  << forceAndEnergySamplingStep*simulationBoxSide << "?\n";
        std::cin >> potentialPrintingStep;
        std::cout << "Do you want to set up a maximum (absolute) value for the prints?\n"
                  << "If yes write it, if no write a negative number. (PLEASE NO ZERO)\n";
        std::cin >> cutoffValue;
        printPotentialsToFileLAMMPS(potentialPrintingStep, cutoffValue);
    }
    else {
        std::cout << "Wrong selection!\n";
    }
}


//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
//************************************************************************//
void IPCsimulation::computeSystemMomentum(double (&pcm)[3]) {
    for (int i: {0, 1, 2})
        pcm[i] = 0.;

    for(IPC ipc: particles) {
        for (int i: {0, 1, 2}) {
            pcm[i] += ipcCenterMass*ipc.ipcCenter.v[i] + firstPatchMass*ipc.firstPatch.v[i] + secndPatchMass*ipc.secndPatch.v[i];
        }
    }
    for(JanusIPC ipc: janusParticles) {
        for (int i: {0, 1, 2}) {
            pcm[i] += ipcCenterMass*ipc.ipcCenter.v[i] + firstPatchMass*ipc.janusPatch.v[i];
        }
    }
}
void IPCsimulation::correctTotalMomentumToZero(double (&pcm)[3], double (&pcmCorrected)[3]) {
    for (int i: {0, 1, 2}) {
        pcmCorrected[i] = 0.;
        pcm[i] /= 3*nIPCs;
    }

    for(IPC ipc: particles) {
        for (int i: {0, 1, 2}) {
            ipc.ipcCenter.v[i]  -= pcm[i];
            ipc.firstPatch.v[i] -= pcm[i];
            ipc.secndPatch.v[i] -= pcm[i];

            pcmCorrected[i] += ipcCenterMass*ipc.ipcCenter.v[i] + firstPatchMass*ipc.firstPatch.v[i] + secndPatchMass*ipc.secndPatch.v[i];
        }
    }
    for(JanusIPC ipc: janusParticles) {
        for (int i: {0, 1, 2}) {
            ipc.ipcCenter.v[i]  -= pcm[i];
            ipc.janusPatch.v[i] -= pcm[i];

            pcmCorrected[i] += ipcCenterMass*ipc.ipcCenter.v[i] + firstPatchMass*ipc.janusPatch.v[i];
        }
    }
}





/*****************************************************************************************/
void IPCsimulation::initializeSystem(const SimulationStage &stage)
{
    isJanusSimulation = stage.janusSimulation;
    printForces = stage.printForces;

    simulationTime = 0;

    initialTemperature = stage.inputStartingTemperature;
    simulationTotalDuration = stage.inputStageTotalDuration;
    printTrajectoryAndCorrelations = stage.inputPrintTrajectoryAndCorrelations;
    readInputFile();

    // if restoring, read state, so we get access to the real number of IPCs
    if(stage.inputRestoringPreviousSimulation) {
        outputFile << "Resuming a previous simulation. ";

        if(isJanusSimulation)
            restorePreviousJanusConfiguration();
        else
            restorePreviousIPCconfiguration();

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
    ratioBetweenTemperatureAndKineticEnergy = 2./(5.*nIPCs-3.);
    ipcRadius = firstPatchEccentricity + firstPatchRadius;  // works for both 2patch and Janus
    interactionRange = 2*ipcRadius;

    // if binary mixture, compute the number of oppositely-charged IPCs
    if(stage.binaryMixturePercentage > 0)
        binaryMixtureComposition = (int) (.01*nIPCs*stage.binaryMixturePercentage);
    else
        binaryMixtureComposition = 0;

    // output the data for future checks
    printInputFileToOutputFile();

    // potential sampling
    outputFile << "Printing potential plots in 'potentials.out'." << std::endl;
    compileForceAndPotentialTables();

    // scale the lenghts to be in a [0.0:1.0] simulation box
    ipcRadius /= simulationBoxSide;
    interactionRange /= simulationBoxSide;
    firstPatchRadius /= simulationBoxSide;
    firstPatchEccentricity /= simulationBoxSide;
    secndPatchRadius /= simulationBoxSide;
    secndPatchEccentricity /= simulationBoxSide;
    dt = simulationTimeStep/simulationBoxSide;
    forceAndEnergySamplingStep /= simulationBoxSide;

    // finish processing data
    squaredInteractionRange = std::pow(interactionRange,2);
    patchDistance = firstPatchEccentricity + secndPatchEccentricity; // the second is zero for Janus
    squaredPatchDistance = patchDistance*patchDistance;
    inversePatchDistance = 1./patchDistance;
    ipcCenterInverseMass = 1./ipcCenterMass;
    firstPatchInverseMass = 1./firstPatchMass;
    halfDtFirstPatchInverseMass = .5*dt*firstPatchInverseMass;
    if(isNotJanusSimulation()) {
        secndPatchInverseMass = 1./secndPatchMass;
        halfDtSecndPatchInverseMass = .5*dt*secndPatchInverseMass;
        // inverse of the I parameter from formulas!
        const double iI = 1./(squaredPatchDistance*ipcCenterInverseMass + secndPatchInverseMass*std::pow(firstPatchEccentricity,2) + firstPatchInverseMass*std::pow(secndPatchEccentricity,2));
        cP11 = 1. - std::pow(secndPatchEccentricity,2)*iI*firstPatchInverseMass;
        cP12 = -firstPatchEccentricity*secndPatchEccentricity*iI*secndPatchInverseMass;
        cP1c = patchDistance*secndPatchEccentricity*iI*ipcCenterInverseMass;
        cP21 = -firstPatchEccentricity*secndPatchEccentricity*iI*firstPatchInverseMass;
        cP22 = 1. - std::pow(firstPatchEccentricity,2)*iI*secndPatchInverseMass;
        cP2c = patchDistance*firstPatchEccentricity*iI*ipcCenterInverseMass;
        alpha_1 = 1. - secndPatchEccentricity*iI*(secndPatchEccentricity*firstPatchInverseMass - firstPatchEccentricity*secndPatchInverseMass);
        alpha_2 = 1. + firstPatchEccentricity*iI*(secndPatchEccentricity*firstPatchInverseMass - firstPatchEccentricity*secndPatchInverseMass);
        alpha_1_firstPatchInverseMass = alpha_1*firstPatchInverseMass;
        alpha_2_secndPatchInverseMass = alpha_2*secndPatchInverseMass;
        alpha_sum = alpha_1_firstPatchInverseMass + alpha_2_secndPatchInverseMass;
        inverseAlpha_sumSquaredPatchDistance = 1./(alpha_sum*squaredPatchDistance);
    }

    // if not restoring, we need to initialize the system here, so that the eccentricities have already been scaled
    if(!stage.inputRestoringPreviousSimulation) {
        outputFile << "Placing " << nIPCs <<  " IPCs on a FCC lattice.\n\n";

        if(isJanusSimulation)
            initializeNewJanusConfiguration();
        else
            initializeNewIPCconfiguration();
    }

    if (binaryMixtureComposition > 0) {
        outputFile << "Binary mixture where " << stage.binaryMixturePercentage << "% of the particles have the opposite charge;\n"
                   << binaryMixtureComposition << " have charge of the opposite sign of that of the other " << nIPCs - binaryMixtureComposition << " particles.\n\n";
    }

    // cell list compilation
    cells.initialize(1., interactionRange, nIPCs);
    outputFile << "Total number of cells: " << cells.getNumberofCells() << std::endl;

    // first computation of forces
    if (isJanusSimulation) {
        cells.compileLists(janusParticles);
        computeFreeJanusForces();
    }
    else {
        cells.compileLists(particles);
        computeFreeForces();
    }

    // check that total momentum is zero
    double pcm [3];
    computeSystemMomentum(pcm);
    outputFile << "P whole system = ( "
               << pcm[0]*simulationBoxSide << ", "
               << pcm[1]*simulationBoxSide << ", "
               << pcm[2]*simulationBoxSide << " )." << std::endl;

    // if not restoring, correct the total momentum to be zero
    double pcmCorrected [3];
    correctTotalMomentumToZero(pcm, pcmCorrected);
    outputFile << "P whole system corrected = ( "
               << pcmCorrected[0]*simulationBoxSide << ", "
               << pcmCorrected[1]*simulationBoxSide << ", "
               << pcmCorrected[2]*simulationBoxSide << " )." << std::endl;

    // first computation of the kinetic energy
    computeSystemEnergy();

    if(stage.inputRestoringPreviousSimulation && initialTemperature > 0) {
        // scale velocities to obtain the desired temperature
        double scalingFactor = std::sqrt(initialTemperature/temperature);
        scaleVelocities(scalingFactor);

        // update energies to include the correction
        computeSystemEnergy();
    }
}

void IPCsimulation::readInputFile() {
    // read input.in file
    std::ifstream inputFile("input.in");
    if(inputFile.fail()) {
        std::cerr << "File input.in could not be opened. Aborting.\n";
        exit(1);
    }
    int N1;
    inputFile >> N1 >> density;
    nIPCs = 4*N1*N1*N1;
    inputFile >> simulationTimeStep >> printingInterval;
    inputFile >> e_BB >> e_Bs1 >> e_Bs2;
    inputFile >> e_s1s1 >> e_s2s2 >> e_s1s2;
    inputFile >> e_min;
    inputFile >> firstPatchEccentricity >> firstPatchRadius;
    inputFile >> secndPatchEccentricity >> secndPatchRadius;
    inputFile >> firstPatchMass >> secndPatchMass >> ipcCenterMass;
    inputFile >> fakeHScoefficient >> fakeHSexponent;
    inputFile >> forceAndEnergySamplingStep >> tollerance;
    inputFile >> isFieldEnabled;
    if(isFieldEnabled) {
        inputFile >> ratioChargeFirstPatchOverIpcCenter >> ratioChargeSecndPatchOverIpcCenter;
        inputFile >> externalFieldIpcCenter[0] >> externalFieldIpcCenter[1] >> externalFieldIpcCenter[2];
        // compute external fields
        for (int i: {0, 1, 2}) {
            externalFieldFirstPatch[i] = ratioChargeFirstPatchOverIpcCenter*externalFieldIpcCenter[i];
            externalFieldSecndPatch[i] = ratioChargeSecndPatchOverIpcCenter*externalFieldIpcCenter[i];
        }
    }
    inputFile.close();
    if(isJanusSimulation) {
        secndPatchRadius = 0.;
        secndPatchEccentricity = 0.;
        secndPatchMass = 0.;
    }

    // patch geometry integrity check
    if ( isNotJanusSimulation() && std::abs( (firstPatchEccentricity+firstPatchRadius)-(secndPatchEccentricity+secndPatchRadius) ) >= 1e-10 ) {
        std::cerr << firstPatchEccentricity << "+" << firstPatchRadius << "=" << firstPatchEccentricity+firstPatchRadius << "-";
        std::cerr << secndPatchEccentricity << "+" << secndPatchRadius << "=" << secndPatchEccentricity+secndPatchRadius << "=\n";
        std::cerr << (firstPatchEccentricity+firstPatchRadius)-(secndPatchEccentricity+secndPatchRadius) << std::endl;
        std::cerr << "eccentricities and radii are not consistent!\n";
        exit(1);
    }
}

void IPCsimulation::printInputFileToOutputFile() {
    outputFile << nIPCs << "\t" << density << "\t" << initialTemperature << "\n";
    outputFile << simulationTimeStep << "\t" << printingInterval << "\t" << simulationTotalDuration << "\n";
    outputFile << e_BB << "\t" << e_Bs1 << "\t" << e_Bs2 << "\n";
    outputFile << e_s1s2 << "\t" << e_s1s1 << "\t" << e_s2s2 << "\n";
    outputFile << e_min << "\n";
    outputFile << firstPatchEccentricity << "\t" << firstPatchRadius << "\n";
    outputFile << secndPatchEccentricity << "\t" << secndPatchRadius << "\n";
    outputFile << firstPatchMass << "\t" << secndPatchMass << "\t" << ipcCenterMass << "\n";
    outputFile << fakeHScoefficient << "\t" << fakeHSexponent << "\n";
    outputFile << forceAndEnergySamplingStep << "\t" << tollerance << "\n";
    outputFile << isFieldEnabled << "\n";
    if(isFieldEnabled) {
        outputFile << ratioChargeFirstPatchOverIpcCenter << "\t" << ratioChargeFirstPatchOverIpcCenter << "\n";
        outputFile << externalFieldIpcCenter[0] << "\t" << externalFieldIpcCenter[1] << "\t" << externalFieldIpcCenter[2] << "\n";
    }

    outputFile << "\n*****************MD simulation in EVN ensemble for CGDH potential.********************\n";
    outputFile << "\nDensity = " << nIPCs << "/" << std::pow(simulationBoxSide,3) << " = ";
    outputFile << nIPCs/std::pow(simulationBoxSide,3) << " = " << density;
    outputFile << "\nSide = " << simulationBoxSide << ", IPC size in reduced units: " << 1./simulationBoxSide << std::endl;
    outputFile << "Total number of sites being simulated: " << (isJanusSimulation? 2*nIPCs : 3*nIPCs) << std::endl;
}

/*****************************************************************************************/


// Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
void IPCsimulation::generateRandomOrientation(double (&a)[3], RandomNumberGenerator & r) {
    double x,y,quad=2.;
    while ( quad > 1. ) {
        x = r.getRandom11();
        y = r.getRandom11();
        quad = x*x + y*y;
    }
    double norm = 2.*sqrt(1.-quad);  a[0]=x*norm;  a[1]=y*norm;  a[2]=1.-2.*quad;
}



//************************************************************************//
double IPCsimulation::computeOmega(double Ra, double Rb, double rab) {
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
double IPCsimulation::computeOmegaRadialDerivative(double Ra, double Rb, double rab) {
    // BKL paper, derivative of formula 18
    if ( rab >= Ra+Rb || rab <= fabs(Ra-Rb) )
        return 0.;
    else {
        const double tempSum = (Ra*Ra-Rb*Rb)/(2.*rab);
        const double tempSumMinus = tempSum - rab/2.;
        const double tempSumPlus = tempSum + rab/2.;
        return (6./rab) * (tempSumMinus*(Ra - tempSumPlus)*(Ra + tempSumPlus) - tempSumPlus*(Rb - tempSumMinus)*(Rb + tempSumMinus) );
    }
}

void IPCsimulation::compileForceAndPotentialTables()
{
    const size_t potentialRangeSamplingSize = size_t( interactionRange/forceAndEnergySamplingStep ) + 1;

    uHS.resize(potentialRangeSamplingSize);
    uBB.resize(potentialRangeSamplingSize);
    uBs1.resize(potentialRangeSamplingSize);
    uBs2.resize(potentialRangeSamplingSize);
    us1s2.resize(potentialRangeSamplingSize);
    us1s1.resize(potentialRangeSamplingSize);
    us2s2.resize(potentialRangeSamplingSize);
    fHS.resize(potentialRangeSamplingSize);
    fBB.resize(potentialRangeSamplingSize);
    fBs1.resize(potentialRangeSamplingSize);
    fBs2.resize(potentialRangeSamplingSize);
    fs1s2.resize(potentialRangeSamplingSize);
    fs1s1.resize(potentialRangeSamplingSize);
    fs2s2.resize(potentialRangeSamplingSize);

    for ( size_t i = 0; i < potentialRangeSamplingSize; ++i)
    {
        const double r = i*forceAndEnergySamplingStep;
        uBB[i]   = (e_BB  /e_min) * computeOmega(ipcRadius, ipcRadius, r);
        uBs1[i]  = (e_Bs1 /e_min) * computeOmega(ipcRadius, firstPatchRadius,  r);
        uBs2[i]  = (e_Bs2 /e_min) * computeOmega(ipcRadius, secndPatchRadius,  r);
        us1s2[i] = (e_s1s2/e_min) * computeOmega(firstPatchRadius,  secndPatchRadius,  r);
        us2s2[i] = (e_s2s2/e_min) * computeOmega(secndPatchRadius,  secndPatchRadius,  r);
        us1s1[i] = (e_s1s1/e_min) * computeOmega(firstPatchRadius,  firstPatchRadius,  r);

        fBB[i]   = (e_BB  /e_min) * computeOmegaRadialDerivative(ipcRadius, ipcRadius, r);
        fBs1[i]  = (e_Bs1 /e_min) * computeOmegaRadialDerivative(ipcRadius, firstPatchRadius,  r);
        fBs2[i]  = (e_Bs2 /e_min) * computeOmegaRadialDerivative(ipcRadius, secndPatchRadius,  r);
        fs1s2[i] = (e_s1s2/e_min) * computeOmegaRadialDerivative(firstPatchRadius,  secndPatchRadius,  r);
        fs2s2[i] = (e_s2s2/e_min) * computeOmegaRadialDerivative(secndPatchRadius,  secndPatchRadius,  r);
        fs1s1[i] = (e_s1s1/e_min) * computeOmegaRadialDerivative(firstPatchRadius,  firstPatchRadius,  r);

        if ( r <= 1.0 )
        {
            // setting up a Fake Hard Sphere Core
            double rm = pow(r, -fakeHSexponent);
            uHS[i]   += fakeHScoefficient*((rm-2.)*rm+1.);
            fHS[i]   += -2.*fakeHSexponent*fakeHScoefficient*(rm-1.)*rm/r;
        }
        // and finally, this division is done here so we don't have to do it during runtime.
        // it comes from the force being Fx = -du/dr dr/dx = -du/dr (x/r)
        const double ir = 1./(r);
        fHS[i]   *= ir;
        fBB[i]   *= ir;
        fBs1[i]  *= ir;
        fBs2[i]  *= ir;
        fs1s2[i] *= ir;
        fs1s1[i] *= ir;
        fs2s2[i] *= ir;
    }
}


void IPCsimulation::computeTrajectoryStep() {
    if (isJanusSimulation) {
        for(JanusIPC &ipc: janusParticles)
            computeVerletHalfStepForJanusIPC(ipc);

        cells.compileLists(janusParticles);
        computeFreeJanusForces();

        for(JanusIPC &ipc: janusParticles)
            finishVerletStepForJanusIPC(ipc);
    }
    else {
        for(IPC &ipc: particles)
            computeVerletHalfStepForIPC(ipc);

        cells.compileLists(particles);
        computeFreeForces();

        for(IPC &ipc: particles)
            finishVerletStepForIPC(ipc);

    }
}


void IPCsimulation::computeSystemEnergy() {
    kineticEnergy = 0.;

    for(IPC ipc: particles) {
        kineticEnergy += firstPatchMass*(std::pow(ipc.firstPatch.v[0],2) + std::pow(ipc.firstPatch.v[1],2) + std::pow(ipc.firstPatch.v[2],2))
           + secndPatchMass*(std::pow(ipc.secndPatch.v[0],2) + std::pow(ipc.secndPatch.v[1],2) + std::pow(ipc.secndPatch.v[2],2))
           + ipcCenterMass*(std::pow(ipc.ipcCenter.v[0],2) + std::pow(ipc.ipcCenter.v[1],2) + std::pow(ipc.ipcCenter.v[2],2));
    }

    for(JanusIPC ipc: janusParticles) {
        kineticEnergy += ipcCenterMass*(std::pow(ipc.ipcCenter.v[0],2) + std::pow(ipc.ipcCenter.v[1],2) + std::pow(ipc.ipcCenter.v[2],2))
                       + firstPatchMass*(std::pow(ipc.janusPatch.v[0],2) + std::pow(ipc.janusPatch.v[1],2) + std::pow(ipc.janusPatch.v[2],2));
    }

    kineticEnergy *= .5*simulationBoxSide*simulationBoxSide;
    totalEnergy = kineticEnergy + potentialEnergy;
    temperature = ratioBetweenTemperatureAndKineticEnergy*kineticEnergy;
}



void IPCsimulation::scaleVelocities(const double scalingFactor) {
    for (IPC &ipc: particles) {
        for (int i: {0, 1, 2}) {
            ipc.ipcCenter.v[i]  *= scalingFactor;
            ipc.firstPatch.v[i] *= scalingFactor;
            ipc.secndPatch.v[i] *= scalingFactor;
        }
    }

    for (JanusIPC &ipc: janusParticles) {
        for (int i: {0, 1, 2}) {
            ipc.ipcCenter.v[i]  *= scalingFactor;
            ipc.janusPatch.v[i] *= scalingFactor;
        }
    }
}


//************************************************************************//
void IPCsimulation::outputSystemTrajectory(std::ofstream & outputTrajectoryFile, const bool printTheForces) {
    outputTrajectoryFile << (isJanusSimulation? 2*nIPCs : 3*nIPCs) << "\n" << simulationBoxSide << "\t" << simulationTime*simulationTimeStep;
    for (IPC ipc: particles) {
        // center
        outputTrajectoryFile << "\n" << ipc.type
                             << "\t" << ipc.ipcCenter.x[0] << "\t" << ipc.ipcCenter.x[1] << "\t" << ipc.ipcCenter.x[2]
                             << "\t" << ipc.ipcCenter.v[0] << "\t" << ipc.ipcCenter.v[1] << "\t" << ipc.ipcCenter.v[2];
        if (printTheForces)     outputTrajectoryFile << "\t" << ipc.ipcCenter.F[0] << "\t" << ipc.ipcCenter.F[1] << "\t" << ipc.ipcCenter.F[2];
        // first patch
        outputTrajectoryFile << "\n" << 'P'
                             << "\t" << ipc.firstPatch.x[0] << "\t" << ipc.firstPatch.x[1] << "\t" << ipc.firstPatch.x[2]
                             << "\t" << ipc.firstPatch.v[0] << "\t" << ipc.firstPatch.v[1] << "\t" << ipc.firstPatch.v[2];
        if (printTheForces)     outputTrajectoryFile << "\t" << ipc.firstPatch.F[0] << "\t" << ipc.firstPatch.F[1] << "\t" << ipc.firstPatch.F[2];
        // second patch
        outputTrajectoryFile << "\n" << 'Q'
                             << "\t" << ipc.secndPatch.x[0] << "\t" << ipc.secndPatch.x[1] << "\t" << ipc.secndPatch.x[2]
                             << "\t" << ipc.secndPatch.v[0] << "\t" << ipc.secndPatch.v[1] << "\t" << ipc.secndPatch.v[2];
        if (printTheForces)     outputTrajectoryFile << "\t" << ipc.secndPatch.F[0] << "\t" << ipc.secndPatch.F[1] << "\t" << ipc.secndPatch.F[2];
    }
    for (JanusIPC ipc: janusParticles) {
        outputTrajectoryFile << "\n" << ipc.type
                             << "\t" << ipc.ipcCenter.x[0] << "\t" << ipc.ipcCenter.x[1] << "\t" << ipc.ipcCenter.x[2]
                             << "\t" << ipc.ipcCenter.v[0] << "\t" << ipc.ipcCenter.v[1] << "\t" << ipc.ipcCenter.v[2];
        if (printTheForces)     outputTrajectoryFile << "\t" << ipc.ipcCenter.F[0] << "\t" << ipc.ipcCenter.F[1] << "\t" << ipc.ipcCenter.F[2];
        outputTrajectoryFile << "\n" << 'P'
                             << "\t" << ipc.janusPatch.x[0] << "\t" << ipc.janusPatch.x[1] << "\t" << ipc.janusPatch.x[2]
                             << "\t" << ipc.janusPatch.v[0] << "\t" << ipc.janusPatch.v[1] << "\t" << ipc.janusPatch.v[2];
        if (printTheForces)     outputTrajectoryFile << "\t" << ipc.janusPatch.F[0] << "\t" << ipc.janusPatch.F[1] << "\t" << ipc.janusPatch.F[2];
    }
    outputTrajectoryFile << std::endl;
}
void IPCsimulation::outputSystemEnergies(std::ofstream & energyTrajectoryFile) {
    energyTrajectoryFile << simulationTime*simulationTimeStep << "\t" << temperature << "\t"
                         << kineticEnergy/nIPCs << "\t" << potentialEnergy/nIPCs << "\t" << totalEnergy/nIPCs << "\t"
                         << std::sqrt(squaredMinimumDistanceBetweenParticles)*simulationBoxSide << std::endl;
}
