#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iomanip>
#include "IPCsimulation.hpp"


//************************************************************************//
IPCsimulation::IPCsimulation(bool restorePreviousSimulation, bool stagingEnabled, std::pair<double,int> stage) {
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
    trajectoryFile.open("siml/trajectory.xyz");
    energyTrajectoryFile.open("siml/evolution.out");

    // initialize system
    initializeSystem(restorePreviousSimulation, stagingEnabled, stage);

    // print starting configuration and initialize output file
    outputFile << "\nPlot evolution.out to check the evolution of the system.\n";

    trajectoryFile << std::scientific << std::setprecision(24);
    energyTrajectoryFile <<std::scientific << std::setprecision(10);
    energyTrajectoryFile << "#t\t\t\tT\t\t\tK\t\t\tU\t\t\tE\t\t\trmin\n";

    outputSystemState(trajectoryFile, energyTrajectoryFile);
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
    int prints = 0;

    // simulation begins
    time(&simulationStartTime);
    while(simulationTime < simulationDurationInIterations) {
        computeTrajectoryStep();
        ++simulationTime;

        if( simulationTime%printingIntervalInIterations == 0) {
            computeSystemEnergy();
            outputSystemState(trajectoryFile, energyTrajectoryFile);
            averageTemperature += temperature;
            averageSquaredTemperature += temperature*temperature;
            averagePotentialEnergy += potentialEnergy;
            ++prints;
        }
    }
    averageTemperature /= prints;
    averagePotentialEnergy /= prints;
    averageSquaredTemperature /= prints;
    double temperatureVariance = std::sqrt(averageSquaredTemperature - std::pow(averageTemperature,2));
    // simulation ends
    time(&simulationEndTime);
    outputFile << "The simulation lasted " << difftime (simulationEndTime, simulationStartTime) << " seconds.\n";

    // output final state
    std::ofstream finalStateFile("startingstate.xyz");
    finalStateFile << std::scientific << std::setprecision(24);
    outputSystemTrajectory(finalStateFile);
    finalStateFile.close();

    // check that total momentum is still zero and print final stuff
    double pcm [3];
    computeSystemMomentum(pcm);
    outputFile << "Residual momentum of the whole system = ( " << pcm[0]*simulationBoxSide << ", " << pcm[1]*simulationBoxSide << ", " << pcm[2]*simulationBoxSide << " ).\n" << std::endl;
    outputFile << "Average kT during the simulation run = " << averageTemperature << std::endl;
    outputFile << "Variance of kT during the simulation run = " << temperatureVariance << std::endl;
    outputFile << "Average potential energy during the simulation run = " << averagePotentialEnergy/nIPCs << std::endl;

    return averageTemperature;
}



//************************************************************************//
void IPCsimulation::computeSystemMomentum(double (&pcm)[3]) {
    for (int i: {0, 1, 2})
        pcm[i] = 0.;

    for(IPC ipc: particles) {
        for (int i: {0, 1, 2}) {
            pcm[i] += ipcCenterMass*ipc.ipcCenter.v[i] + firstPatchMass*ipc.firstPatch.v[i] + secndPatchMass*ipc.secndPatch.v[i];
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
}





/*****************************************************************************************/
void IPCsimulation::initializeSystem(bool restoreprevious, bool stagingEnabled, const std::pair<double,int> & stage)
{
    simulationTime = 0;

    int N1;
    // read input.in file
    std::ifstream inputFile("input.in");
    if(inputFile.fail()) {
        std::cerr << "File input.in could not be opened. Aborting.";
        exit(1);
    }
    inputFile >> N1 >> density >> desiredTemperature;
    nIPCs = 4*N1*N1*N1;
    inputFile >> simulationTimeStep >> printingInterval >> simulationTotalDuration;
    if (stagingEnabled) {
        desiredTemperature = stage.first;
        simulationTotalDuration = stage.second;
    }
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

    // patch geometry integrity check
    if ( abs( (firstPatchEccentricity+firstPatchRadius)-(secndPatchEccentricity+secndPatchRadius) ) >= 1e-10 ) {
        std::cerr << firstPatchEccentricity << "+" << firstPatchRadius << "=" << firstPatchEccentricity+firstPatchRadius << "-";
        std::cerr << secndPatchEccentricity << "+" << secndPatchRadius << "=" << secndPatchEccentricity+secndPatchRadius << "=\n";
        std::cerr << (firstPatchEccentricity+firstPatchRadius)-(secndPatchEccentricity+secndPatchRadius) << std::endl;
        std::cerr << "eccentricities and radii are not consistent!\n";
        exit(1);
    }

    // if restoring, read state, so we get access to the real number of IPCs
    if(restoreprevious) {
        outputFile << "Resuming a previous simulation.\n";
        restorePreviousConfiguration();
        outputFile << "Read " << nIPCs <<  " particles positions and velocities from file.\n\n";
    }

    // process data
    simulationBoxSide = std::cbrt(nIPCs/density);
    ratioBetweenTemperatureAndKineticEnergy = 2./(5.*nIPCs-3.);
    ipcRadius = firstPatchEccentricity + firstPatchRadius;
    interactionRange = 2*ipcRadius;

    // output the data for future checks
    outputFile << nIPCs << "\t" << density << "\t" << desiredTemperature << "\n";
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
    outputFile << "Total number of sites being simulated: " << 3*nIPCs << std::endl;

    // potential sampling
    outputFile << "Printing potential plots in 'potentials.out'." << std::endl;
    make_table(true);

    // scale the lenghts to be in a [0.0:1.0] simulation box
    ipcRadius /= simulationBoxSide;
    interactionRange /= simulationBoxSide;
    firstPatchRadius /= simulationBoxSide;
    firstPatchEccentricity /= simulationBoxSide;
    secndPatchRadius /= simulationBoxSide;
    secndPatchEccentricity /= simulationBoxSide;
    dt = simulationTimeStep/simulationBoxSide;

    // finish processing data
    forceAndEnergySamplingStep /= simulationBoxSide;
    squaredInteractionRange = std::pow(interactionRange,2);
    patchDistance = firstPatchEccentricity + secndPatchEccentricity;
    squaredPatchDistance = patchDistance*patchDistance;
    firstPatchInverseMass = 1./firstPatchMass;
    secndPatchInverseMass = 1./secndPatchMass;
    ipcCenterInverseMass = 1./ipcCenterMass;
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
    alpha_sum = alpha_1 + alpha_2;

    // if not restoring, we need to initialize the system here, so that the eccentricities have already been scaled
    if(!restoreprevious) {
        outputFile << "Starting a new simulation.\n";
        outputFile << "Placing " << nIPCs <<  " IPCs on a FCC lattice.\n\n";
        initializeNewConfiguration(N1);
    }

    // cell list compilation
    cells.initialize(1., interactionRange, nIPCs);
    outputFile << "Total number of cells: " << cells.getNumberofCells() << std::endl;
    cells.compileLists(particles);

    // first computation of forces
    computeFreeForces();

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

    if(restoreprevious && desiredTemperature > 0) {
        // scale velocities to obtain the desired temperature
        double scalingFactor = std::sqrt(desiredTemperature/temperature);
        scaleVelocities(scalingFactor);

        // update energies to include the correction
        computeSystemEnergy();
    }
}





/*****************************************************************************************/


// Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
void IPCsimulation::ranor(double (&a)[3], RandomNumberGenerator & r) {
    double x,y,quad=2.;
    while ( quad > 1. ) {
        x = r.getRandom11();
        y = r.getRandom11();
        quad = x*x + y*y;
    }
    double norm = 2.*sqrt(1.-quad);  a[0]=x*norm;  a[1]=y*norm;  a[2]=1.-2.*quad;
}



//************************************************************************//
double IPCsimulation::omega(double Ra, double Rb, double rab) {
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
double IPCsimulation::d_dr_omega(double Ra, double Rb, double rab) {
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

void IPCsimulation::make_table(bool printPotentials)
{
    const size_t potentialRangeSamplingSize = size_t( interactionRange/forceAndEnergySamplingStep ) + 1;

    uBB.resize(potentialRangeSamplingSize);
    uBs1.resize(potentialRangeSamplingSize);
    uBs2.resize(potentialRangeSamplingSize);
    us1s2.resize(potentialRangeSamplingSize);
    us1s1.resize(potentialRangeSamplingSize);
    us2s2.resize(potentialRangeSamplingSize);
    fBB.resize(potentialRangeSamplingSize);
    fBs1.resize(potentialRangeSamplingSize);
    fBs2.resize(potentialRangeSamplingSize);
    fs1s2.resize(potentialRangeSamplingSize);
    fs1s1.resize(potentialRangeSamplingSize);
    fs2s2.resize(potentialRangeSamplingSize);

    std::ofstream POT_OUTPUT;
    int potOutputPrintCount = 1;
    if (printPotentials) {
        POT_OUTPUT.open("siml/potentials.out");
        POT_OUTPUT << std::scientific << std::setprecision(6);
        POT_OUTPUT << "#r\t\t\tpotBB\t\t\tpotBs1\t\t\tpotBs2\t\t\tpots1s2\t\t\tpots2s2\t\t\tpots1s1";
        POT_OUTPUT <<  "\t\t\tforBB\t\t\tforBs1\t\t\tforBs2\t\t\tfors1s2\t\t\tfors2s2\t\t\tfors1s1\n";
    }

    for ( size_t i = 0; i < potentialRangeSamplingSize; ++i)
    {
        double r = i*forceAndEnergySamplingStep;
        uBB[i]   = (e_BB  /e_min) * omega(ipcRadius, ipcRadius, r);
        uBs1[i]  = (e_Bs1 /e_min) * omega(ipcRadius, firstPatchRadius,  r);
        uBs2[i]  = (e_Bs2 /e_min) * omega(ipcRadius, secndPatchRadius,  r);
        us1s2[i] = (e_s1s2/e_min) * omega(firstPatchRadius,  secndPatchRadius,  r);
        us2s2[i] = (e_s2s2/e_min) * omega(secndPatchRadius,  secndPatchRadius,  r);
        us1s1[i] = (e_s1s1/e_min) * omega(firstPatchRadius,  firstPatchRadius,  r);

        fBB[i]   = (e_BB  /e_min) * d_dr_omega(ipcRadius, ipcRadius, r);
        fBs1[i]  = (e_Bs1 /e_min) * d_dr_omega(ipcRadius, firstPatchRadius,  r);
        fBs2[i]  = (e_Bs2 /e_min) * d_dr_omega(ipcRadius, secndPatchRadius,  r);
        fs1s2[i] = (e_s1s2/e_min) * d_dr_omega(firstPatchRadius,  secndPatchRadius,  r);
        fs2s2[i] = (e_s2s2/e_min) * d_dr_omega(secndPatchRadius,  secndPatchRadius,  r);
        fs1s1[i] = (e_s1s1/e_min) * d_dr_omega(firstPatchRadius,  firstPatchRadius,  r);

        if ( r <= 1.0 )
        {
            // setting up a Fake Hard Sphere Core
            double rm = pow(r, -fakeHSexponent);
            uBB[i]   += fakeHScoefficient*((rm-2.)*rm+1.);
            fBB[i]   += -2.*fakeHSexponent*fakeHScoefficient*(rm-1.)*rm/r;
        }
        if ( printPotentials && int( (1000.*i)/potentialRangeSamplingSize ) == potOutputPrintCount )
        {
            potOutputPrintCount++;
            POT_OUTPUT << r << "\t" << uBB[i] << "\t" << uBs1[i] << "\t" << uBs2[i] << "\t" << us1s2[i] << "\t" << us2s2[i] << "\t" << us1s1[i] << "\t";
            POT_OUTPUT << r << "\t" << fBB[i] << "\t" << fBs1[i] << "\t" << fBs2[i] << "\t" << fs1s2[i] << "\t" << fs2s2[i] << "\t" << fs1s1[i] << "\n";
        }
        // this division is done here to save a division during runtime;
        // it's only done now not to be seen in the plots
        const double x = 1./(r);
        fBB[i]   *= x;
        fBs1[i]  *= x;
        fBs2[i]  *= x;
        fs1s2[i] *= x;
        fs1s1[i] *= x;
        fs2s2[i] *= x;
    }
    POT_OUTPUT.close();
}


void IPCsimulation::computeTrajectoryStep() {
    computeVerletHalfStep();
    cells.compileLists(particles);
    computeFreeForces();
    finishVerletStep();
}

void IPCsimulation::computeVerletHalfStep() {
    for(IPC &ipc: particles) {
        computeVerletHalfStepForIPC(ipc);
    }
}

void IPCsimulation::computeVerletHalfStepForIPC(IPC & ipc) {
    double x1[3], x2[3];
    double dxNew[3];
    for (int i: {0, 1, 2}) {
        // compute the half step velocities from the effective forces of the last step
        ipc.firstPatch.v[i] += ipc.eFp1[i]*(.5*dt*firstPatchInverseMass);
        ipc.secndPatch.v[i] += ipc.eFp2[i]*(.5*dt*secndPatchInverseMass);

        // compute the new positions from the half step velocities
        x1[i] = ipc.firstPatch.x[i] + ipc.firstPatch.v[i]*dt;
        absolutePBC(x1[i]);
        x2[i] = ipc.secndPatch.x[i] + ipc.secndPatch.v[i]*dt;
        absolutePBC(x2[i]);

        // compute the separation between the two patches
        dxNew[i] = x1[i] - x2[i];
        relativePBC(dxNew[i]);
    }
    // compute the (squared) violation of the constraint
    double diff = (dxNew[0]*dxNew[0] + dxNew[1]*dxNew[1] + dxNew[2]*dxNew[2]) - squaredPatchDistance;

    // correct the positions and the velocities until the violation is less than the tollerance
    while( std::fabs(diff) > tollerance*squaredPatchDistance )
    {
        double dxOld[3], DX[3];
        for (int i: {0, 1, 2}) {
            dxOld[i] = ipc.firstPatch.x[i] - ipc.secndPatch.x[i];
            relativePBC(dxOld[i]);
        }
        double g = diff/( 2.*(dxOld[0]*dxNew[0] + dxOld[1]*dxNew[1] + dxOld[2]*dxNew[2]) * alpha_sum*dt );

        for (int i: {0, 1, 2}) {
            DX[i] = g*dxOld[i];

            ipc.firstPatch.v[i] -= alpha_1*DX[i];
            ipc.secndPatch.v[i] += alpha_2*DX[i];

            DX[i] *= dt;

            x1[i] -= DX[i]*alpha_1;
            absolutePBC(x1[i]);
            x2[i] += DX[i]*alpha_2;
            absolutePBC(x2[i]);

            dxNew[i] = x1[i] - x2[i];
            relativePBC(dxNew[i]);
        }

        diff = (dxNew[0]*dxNew[0] + dxNew[1]*dxNew[1] + dxNew[2]*dxNew[2]) - squaredPatchDistance;
    }

    for (int i: {0, 1, 2}) {
        ipc.firstPatch.x[i] = x1[i];
        ipc.secndPatch.x[i] = x2[i];
        ipc.ipcCenter.x[i] = x2[i] + dxNew[i]*secndPatchEccentricity/patchDistance;
        absolutePBC(ipc.ipcCenter.x[i]);
    }
}

void IPCsimulation::finishVerletStep() {
    for(IPC &ipc: particles) {
        finishVerletStepForIPC(ipc);
    }
}

void IPCsimulation::finishVerletStepForIPC(IPC & ipc) {
    double v1[3], v2[3], dx[3], dv[3];
    for (int i: {0, 1, 2}) {
        if (std::fabs(ipc.eFp1[i]) > 1e5) {
            return;
        }
        // compute the the final velocities from the new effective forces
        v1[i] = ipc.firstPatch.v[i] + ipc.eFp1[i]*(.5*dt*firstPatchInverseMass);
        v2[i] = ipc.secndPatch.v[i] + ipc.eFp2[i]*(.5*dt*secndPatchInverseMass);

        // compute the patch-patch distance
        dx[i] = ipc.firstPatch.x[i] - ipc.secndPatch.x[i];
        relativePBC(dx[i]);
        dv[i] = v1[i] - v2[i];
    }
    // check how much the constraints are being violated
    double k = (dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2])/(alpha_sum*squaredPatchDistance);
    while( std::fabs(k) > tollerance ) {
        // compute and apply corrections
        double DX[3];
        for (int i: {0, 1, 2}) {
            DX[i] = k*dx[i];
            v1[i] -= DX[i]*alpha_1;
            v2[i] += DX[i]*alpha_2;
            dv[i] = v1[i] - v2[i];
        }
        // recompute the violation of the constraints
        k = (dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2])/(alpha_sum*squaredPatchDistance);
    }

    for (int i: {0, 1, 2}) {
        ipc.firstPatch.v[i] = v1[i];
        ipc.secndPatch.v[i] = v2[i];
        ipc.ipcCenter.v[i] = (v1[i]*secndPatchEccentricity + v2[i]*firstPatchEccentricity)/patchDistance;
    }
}



void IPCsimulation::computeSystemEnergy() {
    kineticEnergy = 0.;
    for(IPC ipc: particles) {
        kineticEnergy += firstPatchMass*(std::pow(ipc.firstPatch.v[0],2) + std::pow(ipc.firstPatch.v[1],2) + std::pow(ipc.firstPatch.v[2],2))
           + secndPatchMass*(std::pow(ipc.secndPatch.v[0],2) + std::pow(ipc.secndPatch.v[1],2) + std::pow(ipc.secndPatch.v[2],2))
           + ipcCenterMass*(std::pow(ipc.ipcCenter.v[0],2) + std::pow(ipc.ipcCenter.v[1],2) + std::pow(ipc.ipcCenter.v[2],2));
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
}


//************************************************************************//
void IPCsimulation::outputSystemTrajectory(std::ofstream & outputTrajectoryFile) {
    outputTrajectoryFile<<3*nIPCs<<"\n"<<simulationTime*simulationTimeStep;
    for (IPC ipc: particles) {
        outputTrajectoryFile << "\n"
                             << ipc.type << "\t" << ipc.ipcCenter.x[0] << "\t" << ipc.ipcCenter.x[1] << "\t" << ipc.ipcCenter.x[2]
                                         << "\t" << ipc.ipcCenter.v[0] << "\t" << ipc.ipcCenter.v[1] << "\t" << ipc.ipcCenter.v[2]
                                         //<< "\t" << ipc.ipcCenter.F[0] << "\t" << ipc.ipcCenter.F[1] << "\t" << ipc.ipcCenter.F[2]
                             << "\n"
                             << 'P'      << "\t" << ipc.firstPatch.x[0] << "\t" << ipc.firstPatch.x[1] << "\t" << ipc.firstPatch.x[2]
                                         << "\t" << ipc.firstPatch.v[0] << "\t" << ipc.firstPatch.v[1] << "\t" << ipc.firstPatch.v[2]
                                         //<< "\t" << ipc.firstPatch.F[0] << "\t" << ipc.firstPatch.F[1] << "\t" << ipc.firstPatch.F[2]
                             << "\n"
                             << 'Q'      << "\t" << ipc.secndPatch.x[0] << "\t" << ipc.secndPatch.x[1] << "\t" << ipc.secndPatch.x[2]
                                         << "\t" << ipc.secndPatch.v[0] << "\t" << ipc.secndPatch.v[1] << "\t" << ipc.secndPatch.v[2];
                                         //<< "\t" << ipc.secndPatch.F[0] << "\t" << ipc.secndPatch.F[1] << "\t" << ipc.secndPatch.F[2];
    }
    outputTrajectoryFile << std::endl;
}
void IPCsimulation::outputSystemState(std::ofstream & outputTrajectoryFile, std::ofstream & energyTrajectoryFile)
{
    outputSystemTrajectory(outputTrajectoryFile);
    energyTrajectoryFile << simulationTime*simulationTimeStep << "\t" << temperature << "\t"
                         << kineticEnergy/nIPCs << "\t" << potentialEnergy/nIPCs << "\t" << totalEnergy/nIPCs << "\t"
                         << std::sqrt(squaredMinimumDistanceBetweenParticles)*simulationBoxSide << std::endl;
}




void IPCsimulation::initializeNewConfiguration(int N1) {
    particles.resize(nIPCs);
    RandomNumberGenerator rand;

    int N2 = N1*N1;
    int N3 = N2*N1;

    // scaling: sqrt(2kT/mPI) comes from boltzmann average of |v_x|
    //double vel_scaling = std::sqrt(2.*desiredTemperature/3.1415)/simulationBoxSide;

    double vel_scaling = std::sqrt(1.6*desiredTemperature)/simulationBoxSide;
    // initialize IPC positions
    for(int i=0;i<N3;i++)
    {
      // FCC is obtained as 4 intersecating SC
        particles[i].number = i;
        particles[i].type = 'C';
        particles[i].ipcCenter.x[0] = (i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i].ipcCenter.x[0]);
        particles[i].ipcCenter.x[1] = ((i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i].ipcCenter.x[1]);
        particles[i].ipcCenter.x[2] = (i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(particles[i].ipcCenter.x[2]);

        particles[i+N3].number = i+N3;
        particles[i+N3].type = 'C';
        particles[i+N3].ipcCenter.x[0] = (.5 + i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3].ipcCenter.x[0]);
        particles[i+N3].ipcCenter.x[1] = (.5 + (i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3].ipcCenter.x[1]);
        particles[i+N3].ipcCenter.x[2] = (i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(particles[i+N3].ipcCenter.x[2]);

        particles[i+N3+N3].number = i+N3+N3;
        particles[i+N3+N3].type = 'C';
        particles[i+N3+N3].ipcCenter.x[0] = (i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3+N3].ipcCenter.x[0]);
        particles[i+N3+N3].ipcCenter.x[1] = (.5 + (i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3+N3].ipcCenter.x[1]);
        particles[i+N3+N3].ipcCenter.x[2] = (.5 + i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(particles[i+N3+N3].ipcCenter.x[2]);

        particles[i+N3+N3+N3].number = i+N3+N3+N3;
        particles[i+N3+N3+N3].type = 'C';
        particles[i+N3+N3+N3].ipcCenter.x[0] = (.5 + i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3+N3+N3].ipcCenter.x[0]);
        particles[i+N3+N3+N3].ipcCenter.x[1] = ((i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(particles[i+N3+N3+N3].ipcCenter.x[1]);
        particles[i+N3+N3+N3].ipcCenter.x[2] = (.5 + i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(particles[i+N3+N3+N3].ipcCenter.x[2]);

        // starting from random but ONLY TRANSLATIONAL speeds, for compatibility with rattle
        for (int j: {0, 1, 2}) {
           particles[i].ipcCenter.v[j]          = rand.getRandom11()*vel_scaling;
           particles[i+N3].ipcCenter.v[j]       = rand.getRandom11()*vel_scaling;
           particles[i+N3+N3].ipcCenter.v[j]    = rand.getRandom11()*vel_scaling;
           particles[i+N3+N3+N3].ipcCenter.v[j] = rand.getRandom11()*vel_scaling;
        }
    }
    // initialize patches positions
    for(IPC &ipc: particles) {
        double ipcAxis[3], ipcOrthogonalAxis[3];
        ranor(ipcAxis,rand);
        ranor(ipcOrthogonalAxis, rand);
        double normOfIpcOrthogonalAxis = std::sqrt( .5*
                (std::pow(ipc.ipcCenter.v[0],2) + std::pow(ipc.ipcCenter.v[1],2) + std::pow(ipc.ipcCenter.v[2],2))
                /
                (std::pow(ipcOrthogonalAxis[0],2) + std::pow(ipcOrthogonalAxis[1],2) + std::pow(ipcOrthogonalAxis[2],2))
                );
        double scalarProductOfTheTwoAxes = ipcOrthogonalAxis[0]*ipcAxis[0] +
                ipcOrthogonalAxis[1]*ipcAxis[1] + ipcOrthogonalAxis[2]*ipcAxis[2];

        for (int i: {0, 1, 2}) {
            ipcOrthogonalAxis[i] *= normOfIpcOrthogonalAxis;
            ipcOrthogonalAxis[i] -= ipcAxis[i]*scalarProductOfTheTwoAxes;

            ipc.firstPatch.x[i] = ipc.ipcCenter.x[i] + ipcAxis[i]*firstPatchEccentricity;
            absolutePBC(ipc.firstPatch.x[i]);
            ipc.secndPatch.x[i] = ipc.ipcCenter.x[i] - ipcAxis[i]*secndPatchEccentricity;
            absolutePBC(ipc.secndPatch.x[i]);

            double temp = rand.getRandom11()*ipcOrthogonalAxis[i]*vel_scaling;
            ipc.firstPatch.v[i] = ipc.ipcCenter.v[i] + temp;
            ipc.secndPatch.v[i] = ipc.ipcCenter.v[i] - temp;
        }
    }
}

void IPCsimulation::restorePreviousConfiguration() {
    char unusedPatchName;
    double unusedTime;
    std::ifstream startingConfigurationFile("startingstate.xyz");
    if(startingConfigurationFile.fail()) {
        std::cerr << "File startingstate.xyz could not be opened. Aborting.";
        exit(1);
    }
    startingConfigurationFile >> nIPCs >> unusedTime;
    nIPCs /= 3;

    particles.resize(nIPCs);
    int counter = 0;

    for (IPC &ipc: particles) {
        ipc.number = counter++;
        startingConfigurationFile >> ipc.type
           >> ipc.ipcCenter.x[0] >> ipc.ipcCenter.x[1] >> ipc.ipcCenter.x[2]
           >> ipc.ipcCenter.v[0] >> ipc.ipcCenter.v[1] >> ipc.ipcCenter.v[2];
        startingConfigurationFile >> unusedPatchName
           >> ipc.firstPatch.x[0] >> ipc.firstPatch.x[1] >> ipc.firstPatch.x[2]
           >> ipc.firstPatch.v[0] >> ipc.firstPatch.v[1] >> ipc.firstPatch.v[2];
        startingConfigurationFile >> unusedPatchName
           >> ipc.secndPatch.x[0] >> ipc.secndPatch.x[1] >> ipc.secndPatch.x[2]
           >> ipc.secndPatch.v[0] >> ipc.secndPatch.v[1] >> ipc.secndPatch.v[2];
    }
    if (counter != nIPCs) {
        std::cerr << "Placed " << counter << " IPCs, expected " << nIPCs << ", quitting.\n";
        exit(1);
    }

    startingConfigurationFile.close();
}



void IPCsimulation::computeFreeForces() {
    // reset all forces
    for(IPC &ipc: particles) {
        for (int i: {0, 1, 2}) {
            ipc.ipcCenter.F[i] = 0.;
            ipc.firstPatch.F[i] = 0.;
            ipc.secndPatch.F[i] = 0.;
        }
    }
    potentialEnergy = 0.0;  squaredMinimumDistanceBetweenParticles = 1.;

    #pragma omp parallel
    {
        loopVariables loopVars(nIPCs);

        #pragma omp for
        for(int m=0; m<cells.getNumberofCells(); ++m)  // loop over all cells
        {
            const std::list<int> & ipcsInCell = cells.getIPCsInCell(m);
            const std::list<int> ipcsInNeighbouringCells = cells.getIPCsInNeighbouringCells(m);
            for(auto ipc = ipcsInCell.cbegin(); ipc != ipcsInCell.cend(); ++ipc) {
                computeInteractionsWithIPCsInTheSameCell(ipc, ipcsInCell, loopVars);
                computeInteractionsWithIPCsInNeighbouringCells(ipc, ipcsInNeighbouringCells, loopVars);
            }
        }
        #pragma omp critical
        {
            for (IPC &ipc: particles) {
                for (int i: {0, 1, 2}) {
                    const size_t j = ipc.number;
                    ipc.ipcCenter.F[i]  += loopVars.ipcCenterF[j][i];
                    ipc.firstPatch.F[i] += loopVars.firstPatchF[j][i];
                    ipc.secndPatch.F[i] += loopVars.secndPatchF[j][i];
                }
            }
            potentialEnergy += loopVars.U;
            squaredMinimumDistanceBetweenParticles += loopVars.minimumSquaredDistance;
        }
//        #pragma omp for
//        for(int n=0; n < nIPCs; ++n) {
//            IPC& ipc = particles[n];
        for(IPC &ipc: particles) {
            for (int i: {0, 1, 2}) {
                if(isFieldEnabled) {
                    ipc.ipcCenter.F[i]  += externalFieldIpcCenter[i];
                    ipc.firstPatch.F[i] += externalFieldFirstPatch[i];
                    ipc.secndPatch.F[i] += externalFieldSecndPatch[i];
                }
                ipc.eFp1[i] = ipc.firstPatch.F[i]*cP11 + ipc.secndPatch.F[i]*cP12 + ipc.ipcCenter.F[i]*cP1c;
                ipc.eFp2[i] = ipc.firstPatch.F[i]*cP21 + ipc.secndPatch.F[i]*cP22 + ipc.ipcCenter.F[i]*cP2c;
            }
        }
    }
}

void IPCsimulation::computeInteractionsWithIPCsInNeighbouringCells(std::list<int>::const_iterator loc, std::list<int> const& ipcsInNeighbouringCells, loopVariables & loopVars) {
    for( auto ext = ipcsInNeighbouringCells.cbegin(); ext != ipcsInNeighbouringCells.cend(); ++ext) {
        computeInteractionsBetweenTwoIPCs(*loc, *ext, loopVars);
    }
}



void IPCsimulation::computeInteractionsWithIPCsInTheSameCell(std::list<int>::const_iterator loc, std::list<int> const& ipcsInCurrentCell, loopVariables &loopVars) {
    // starts from loc+1 which is like summing over i > j inside the cell
    for(std::list<int>::const_iterator ins = std::next(loc); ins != ipcsInCurrentCell.cend(); ++ins) {
        computeInteractionsBetweenTwoIPCs(*loc, *ins, loopVars);
    }
}

void IPCsimulation::computeInteractionsBetweenTwoIPCs(const int firstIPC, const int secndIPC, loopVariables &loopVars) {

    IPC const& first = particles[firstIPC];
    IPC const& secnd = particles[secndIPC];

    // compute center-center distance
    double centerCenterSeparation[3];
    for (int i: {0, 1, 2}) {
        centerCenterSeparation[i] = first.ipcCenter.x[i] - secnd.ipcCenter.x[i];
        relativePBC(centerCenterSeparation[i]);
    }
    double centerCenterSeparationModulus = centerCenterSeparation[0]*centerCenterSeparation[0]
                                         + centerCenterSeparation[1]*centerCenterSeparation[1]
                                         + centerCenterSeparation[2]*centerCenterSeparation[2];

    if (centerCenterSeparationModulus < loopVars.minimumSquaredDistance)
        loopVars.minimumSquaredDistance = centerCenterSeparationModulus;

    // if the CENTERS are too far, no interactions, skip this couple of IPCs
    if (centerCenterSeparationModulus >= squaredInteractionRange)
        return;

    // we are inside the interaction range; compute the interaction between centers
    centerCenterSeparationModulus = std::sqrt(centerCenterSeparationModulus);
    const size_t centerCenterDistance = size_t( centerCenterSeparationModulus/forceAndEnergySamplingStep );
    loopVars.U += uBB[centerCenterDistance];
    for (int i: {0, 1, 2}) {
        const double modulus = fBB[centerCenterDistance]*centerCenterSeparation[i];
        loopVars.ipcCenterF[firstIPC][i] -= modulus;
        loopVars.ipcCenterF[secndIPC][i] += modulus;
    }

    // compute all the other site-site separations
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

    // all the others
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
            loopVars.U += uBs1[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fBs1[dist]*siteSiteSeparation[0][i];
                loopVars.ipcCenterF[firstIPC][i] -= modulus;
                loopVars.firstPatchF[secndIPC][i] += modulus;
            }
        } else if (j == 1) { // center - patch2
            loopVars.U += uBs2[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fBs2[dist]*siteSiteSeparation[1][i];
                loopVars.ipcCenterF[firstIPC][i] -= modulus;
                loopVars.secndPatchF[secndIPC][i] += modulus;
            }
        } else if (j == 2) { // patch1 - center
            loopVars.U += uBs1[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fBs1[dist]*siteSiteSeparation[2][i];
                loopVars.firstPatchF[firstIPC][i] -= modulus;
                loopVars.ipcCenterF[secndIPC][i] += modulus;
            }
        } else if (j == 3) { // patch1 - patch1
            loopVars.U += us1s1[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fs1s1[dist]*siteSiteSeparation[3][i];
                loopVars.firstPatchF[firstIPC][i] -= modulus;
                loopVars.firstPatchF[secndIPC][i] += modulus;
            }
        } else if (j == 4) { // patch1 - patch2
            loopVars.U += us1s2[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fs1s2[dist]*siteSiteSeparation[4][i];
                loopVars.firstPatchF[firstIPC][i] -= modulus;
                loopVars.secndPatchF[secndIPC][i] += modulus;
            }
        } else if (j == 5) { // patch2 - center
            loopVars.U += uBs2[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fBs2[dist]*siteSiteSeparation[5][i];
                loopVars.secndPatchF[firstIPC][i] -= modulus;
                loopVars.ipcCenterF[secndIPC][i] += modulus;
            }
        } else if (j == 6) { // patch2 - patch1
            loopVars.U += us1s2[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fs1s2[dist]*siteSiteSeparation[6][i];
                loopVars.secndPatchF[firstIPC][i] -= modulus;
                loopVars.firstPatchF[secndIPC][i] += modulus;
            }
        } else if (j == 7) { // patch2 - patch2
            loopVars.U += us2s2[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fs2s2[dist]*siteSiteSeparation[7][i];
                loopVars.secndPatchF[firstIPC][i] -= modulus;
                loopVars.secndPatchF[secndIPC][i] += modulus;
            }
        }
    }
}
