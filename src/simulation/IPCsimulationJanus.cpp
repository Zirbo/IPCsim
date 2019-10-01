#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iomanip>
#include "IPCsimulation.hpp"





/*****************************************************************************************/
void JanusIPCsimulation::initializeSystem(bool restoreprevious, bool stagingEnabled, const std::pair<double,int> & stage)
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
    inputFile >> e_BB >> e_Bs >> e_ss;
    inputFile >> e_min;
    inputFile >> patchCenterDistance >> patchRadius;
    inputFile >> patchMass >> centerMass;
    inputFile >> fakeHScoefficient >> fakeHSexponent;
    inputFile >> forceAndEnergySamplingStep >> tollerance;
    inputFile >> isFieldEnabled;
    if(isFieldEnabled) {
        inputFile >> ratioChargePatchOverIpcCenter;
        inputFile >> externalFieldIpcCenter[0] >> externalFieldIpcCenter[1] >> externalFieldIpcCenter[2];
        // compute external fields
        for (int i: {0, 1, 2}) {
            externalFieldPatch[i] = ratioChargePatchOverIpcCenter*externalFieldIpcCenter[i];
        }
    }
    inputFile.close();

    // patch geometry integrity check
    if ( patchCenterDistance >= .5 ) {
        std::cerr << "The patch-center distance is " << patchCenterDistance
                  << ", bigger than the radius of the hard core!" << std::endl;
        exit(1);
    }

    // if restoring, read state, so we get access to the real number of IPCs
    if(restoreprevious) {
        outputFile << "Resuming a previous simulation.\n";
        restorePreviousConfiguration();
        outputFile << "Read " << nIPCs <<  " particles positions and velocities from file.\n\n";
        // we read nIPCs and simulationBoxSide from the starting configuration, so we can compute the density from them
        density = double(nIPCs)/std::pow(simulationBoxSide, 3);
    }
    else {
        // we read nIPCs and density from the input file, so we need to compute the simulationBoxSide from them
        simulationBoxSide = std::cbrt(nIPCs/density);
    }

    // process data
    simulationBoxSide = std::cbrt(nIPCs/density);
    ratioBetweenTemperatureAndKineticEnergy = 2./(5.*nIPCs-3.);
    ipcRadius = patchCenterDistance + patchRadius;
    interactionRange = 2*ipcRadius;

    // output the data for future checks
    outputFile << nIPCs << "\t" << density << "\t" << desiredTemperature << "\n";
    outputFile << simulationTimeStep << "\t" << printingInterval << "\t" << simulationTotalDuration << "\n";
    outputFile << e_BB << "\t" << e_Bs << "\t" << e_ss << "\n";
    outputFile << e_min << "\n";
    outputFile << patchCenterDistance << "\t" << patchRadius << "\n";
    outputFile << patchMass << "\t" << centerMass << "\n";
    outputFile << fakeHScoefficient << "\t" << fakeHSexponent << "\n";
    outputFile << forceAndEnergySamplingStep << "\t" << tollerance << "\n";
    outputFile << isFieldEnabled << "\n";
    if(isFieldEnabled) {
        outputFile << ratioChargePatchOverIpcCenter << "\t" << ratioChargePatchOverIpcCenter << "\n";
        outputFile << externalFieldIpcCenter[0] << "\t" << externalFieldIpcCenter[1] << "\t" << externalFieldIpcCenter[2] << "\n";
    }

    outputFile << "\n*****************MD simulation in EVN ensemble for CGDH potential.********************\n";
    outputFile << "\nDensity = " << nIPCs << "/" << std::pow(simulationBoxSide,3) << " = ";
    outputFile << nIPCs/std::pow(simulationBoxSide,3) << " = " << density;
    outputFile << "\nSide = " << simulationBoxSide << ", IPC size in reduced units: " << 1./simulationBoxSide << std::endl;
    outputFile << "Total number of sites being simulated: " << 2*nIPCs << std::endl;

    // potential sampling
    outputFile << "Printing potential plots in 'potentials.out'." << std::endl;
    make_table(true);

    // scale the lenghts to be in a [0.0:1.0] simulation box
    ipcRadius /= simulationBoxSide;
    interactionRange /= simulationBoxSide;
    patchRadius /= simulationBoxSide;
    patchCenterDistance /= simulationBoxSide;
    dt = simulationTimeStep/simulationBoxSide;

    // finish processing data
    forceAndEnergySamplingStep /= simulationBoxSide;
    squaredInteractionRange = std::pow(interactionRange,2);
    patchInverseMass = 1./patchMass;
    centerInverseMass = 1./centerMass;

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



void JanusIPCsimulation::make_table(bool printPotentials)
{
    const size_t potentialRangeSamplingSize = size_t( interactionRange/forceAndEnergySamplingStep ) + 1;

    uBB.resize(potentialRangeSamplingSize);
    uBs.resize(potentialRangeSamplingSize);
    uss.resize(potentialRangeSamplingSize);
    fBB.resize(potentialRangeSamplingSize);
    fBs.resize(potentialRangeSamplingSize);
    fss.resize(potentialRangeSamplingSize);

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
        uBB[i]  = (e_BB /e_min) * omega(ipcRadius, ipcRadius, r);
        uBs[i]  = (e_Bs /e_min) * omega(ipcRadius, patchRadius, r);
        uss[i]  = (e_ss /e_min) * omega(patchRadius, patchRadius, r);

        fBB[i]  = (e_BB /e_min) * d_dr_omega(ipcRadius, ipcRadius, r);
        fBs[i]  = (e_Bs /e_min) * d_dr_omega(ipcRadius, patchRadius, r);
        fss[i]  = (e_ss /e_min) * d_dr_omega(patchRadius, patchRadius, r);

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
            POT_OUTPUT << r << "\t" << uBB[i] << "\t" << uBs[i] << "\t" << uss[i] << "\t";
            POT_OUTPUT << r << "\t" << fBB[i] << "\t" << fBs[i] << "\t" << fss[i] << "\n";
        }
        // this division is done here to save a division during runtime;
        // it's only done now not to be seen in the plots
        const double x = 1./(r);
        fBB[i] *= x;
        fBs[i] *= x;;
        fss[i] *= x;
    }
    POT_OUTPUT.close();
}

void IPCsimulation::computeVerletHalfStepForIPCJanus(JanusIPC & ipc) {
    double xc[3], xp[3];
    double dxNew[3];
    for (int i: {0, 1, 2}) {
        // compute the half step velocities from the effective forces of the last step
        ipc.ipcCenter.v[i] += ipc.ipcCenter.F[i]*(.5*dt*centerInverseMass);
        ipc.janusPatch.v[i] += ipc.janusPatch.F[i]*(.5*dt*patchInverseMass);

        // compute the new positions from the half step velocities
        xc[i] = ipc.ipcCenter.x[i] + ipc.ipcCenter.v[i]*dt;
        absolutePBC(xc[i]);
        xp[i] = ipc.janusPatch.x[i] + ipc.janusPatch.v[i]*dt;
        absolutePBC(xp[i]);

        // compute the separation between the two patches
        dxNew[i] = xc[i] - xp[i];
        relativePBC(dxNew[i]);
    }
    // compute the (squared) violation of the constraint
    double diff = (dxNew[0]*dxNew[0] + dxNew[1]*dxNew[1] + dxNew[2]*dxNew[2]) - squaredPatchCenterDistance;

    // correct the positions and the velocities until the violation is less than the tollerance
    while( std::fabs(diff) > tollerance*squaredPatchCenterDistance )
    {
        double dxOld[3], DX[3];
        for (int i: {0, 1, 2}) {
            dxOld[i] = ipc.ipcCenter.x[i] - ipc.janusPatch.x[i];
            relativePBC(dxOld[i]);
        }
        double g = diff/( 2.*(dxOld[0]*dxNew[0] + dxOld[1]*dxNew[1] + dxOld[2]*dxNew[2]) * dt );

        for (int i: {0, 1, 2}) {
            DX[i] = g*dxOld[i];

            ipc.ipcCenter.v[i] += DX[i];
            ipc.janusPatch.v[i] -= DX[i];

            DX[i] *= dt;

            xc[i] -= DX[i];
            absolutePBC(xc[i]);
            xp[i] += DX[i];
            absolutePBC(xp[i]);

            dxNew[i] = xc[i] - xp[i];
            relativePBC(dxNew[i]);
        }

        diff = (dxNew[0]*dxNew[0] + dxNew[1]*dxNew[1] + dxNew[2]*dxNew[2]) - squaredPatchCenterDistance;
    }

    for (int i: {0, 1, 2}) {
        ipc.ipcCenter.x[i] = xc[i];
        ipc.janusPatch.x[i] = xp[i];
    }
}


void IPCsimulation::finishVerletStepForIPCJanus(JanusIPC & ipc) {
    double vc[3], vp[3], dx[3], dv[3];
    for (int i: {0, 1, 2}) {
        // compute the the final velocities from the new effective forces
        vc[i] = ipc.ipcCenter.v[i] + ipc.ipcCenter.F[i]*(.5*dt*centerInverseMass);
        vp[i] = ipc.janusPatch.v[i] + ipc.janusPatch.F[i]*(.5*dt*patchInverseMass);

        // compute the patch-patch distance
        dx[i] = ipc.ipcCenter.x[i] - ipc.janusPatch.x[i];
        relativePBC(dx[i]);
        dv[i] = vc[i] - vp[i];
    }
    // check how much the constraints are being violated
    double k = (dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2]) / squaredPatchCenterDistance;
    while( std::fabs(k) > tollerance ) {
        // compute and apply corrections
        double DX[3];
        for (int i: {0, 1, 2}) {
            DX[i] = k*dx[i];
            vc[i] -= DX[i];
            vp[i] += DX[i];
            dv[i] = vc[i] - vp[i];
        }
        // recompute the violation of the constraints
        k = (dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2]) / squaredPatchCenterDistance;
    }

    for (int i: {0, 1, 2}) {
        ipc.ipcCenter.v[i] = vc[i];
        ipc.janusPatch.v[i] = vp[i];
    }
}




void JanusIPCsimulation::initializeNewConfiguration(int N1) {
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
    for(JanusIPC &ipc: particles) {
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

            ipc.janusPatch.x[i] = ipc.ipcCenter.x[i] + ipcAxis[i]*patchCenterDistance;
            absolutePBC(ipc.janusPatch.x[i]);

            double temp = rand.getRandom11()*ipcOrthogonalAxis[i]*vel_scaling;
            ipc.janusPatch.v[i] = ipc.ipcCenter.v[i] + temp;
        }
    }
}

void JanusIPCsimulation::restorePreviousConfiguration() {
    char unusedPatchName;
    double unusedTime;
    std::ifstream startingConfigurationFile("startingstate.xyz");
    if(startingConfigurationFile.fail()) {
        std::cerr << "File startingstate.xyz could not be opened. Aborting.";
        exit(1);
    }
    startingConfigurationFile >> nIPCs >> simulationBoxSide >> unusedTime;
    nIPCs /= 3;

    particles.resize(nIPCs);
    int counter = 0;

    for (JanusIPC &ipc: particles) {
        ipc.number = counter++;
        startingConfigurationFile >> ipc.type
           >> ipc.ipcCenter.x[0] >> ipc.ipcCenter.x[1] >> ipc.ipcCenter.x[2]
           >> ipc.ipcCenter.v[0] >> ipc.ipcCenter.v[1] >> ipc.ipcCenter.v[2];
        startingConfigurationFile >> unusedPatchName
           >> ipc.janusPatch.x[0] >> ipc.janusPatch.x[1] >> ipc.janusPatch.x[2]
           >> ipc.janusPatch.v[0] >> ipc.janusPatch.v[1] >> ipc.janusPatch.v[2];
    }
    if (counter != nIPCs) {
        std::cerr << "Placed " << counter << " IPCs, expected " << nIPCs << ", quitting.\n";
        exit(1);
    }

    startingConfigurationFile.close();
}



void JanusIPCsimulation::computeFreeForces() {
    // reset all forces
    for(JanusIPC &ipc: particles) {
        for (int i: {0, 1, 2}) {
            ipc.ipcCenter.F[i] = 0.;
            ipc.janusPatch.F[i] = 0.;
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
            for (JanusIPC &ipc: particles) {
                for (int i: {0, 1, 2}) {
                    const size_t j = ipc.number;
                    ipc.ipcCenter.F[i]  += loopVars.ipcCenterF[j][i];
                    ipc.janusPatch.F[i] += loopVars.janusPatchF[j][i];
                }
            }
            potentialEnergy += loopVars.U;
            squaredMinimumDistanceBetweenParticles += loopVars.minimumSquaredDistance;
        }
//        #pragma omp for
//        for(int n=0; n < nIPCs; ++n) {
//            IPC& ipc = particles[n];
        for(JanusIPC &ipc: particles) {
            for (int i: {0, 1, 2}) {
                if(isFieldEnabled) {
                    ipc.ipcCenter.F[i] += externalFieldIpcCenter[i];
                    ipc.janusPatch.F[i] += externalFieldPatch[i];
                }
            }
        }
    }
}

void IPCsimulation::computeInteractionsBetweenTwoIPCsJanus(const int firstIPC, const int secndIPC, loopVariables &loopVars) {

    JanusIPC const& first = particles[firstIPC];
    JanusIPC const& secnd = particles[secndIPC];

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
    double siteSiteSeparation[3][3];
    for (int i: {0, 1, 2}) {
        siteSiteSeparation[0][i] = first.ipcCenter.x[i] - secnd.janusPatch.x[i];
        siteSiteSeparation[1][i] = first.janusPatch.x[i] - secnd.ipcCenter.x[i];
        siteSiteSeparation[2][i] = first.janusPatch.x[i] - secnd.janusPatch.x[i];
        for (int j = 0; j < 3; ++j)
            relativePBC(siteSiteSeparation[j][i]);
    }

    // all the others
    for (int j = 0; j < 3; ++j) {
        double siteSiteSeparationModulus = siteSiteSeparation[j][0]*siteSiteSeparation[j][0]
                                         + siteSiteSeparation[j][1]*siteSiteSeparation[j][1]
                                         + siteSiteSeparation[j][2]*siteSiteSeparation[j][2];

        // if we are too far, no interaction, skip to the next site-site pair
        if (siteSiteSeparationModulus >= squaredInteractionRange)
            continue;

        siteSiteSeparationModulus = std::sqrt(siteSiteSeparationModulus);
        const size_t dist = size_t( siteSiteSeparationModulus/forceAndEnergySamplingStep );
        if (j == 0) { // center - patch
            loopVars.U += uBs[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fBs[dist]*siteSiteSeparation[0][i];
                loopVars.ipcCenterF[firstIPC][i] -= modulus;
                loopVars.janusPatchF[secndIPC][i] += modulus;
            }
        } else if (j == 1) { // patch - center
            loopVars.U += uBs[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fBs[dist]*siteSiteSeparation[1][i];
                loopVars.janusPatchF[firstIPC][i] -= modulus;
                loopVars.ipcCenterF[secndIPC][i] += modulus;
            }
        } else if (j == 2) { // patch - patch
            loopVars.U += uss[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fss[dist]*siteSiteSeparation[2][i];
                loopVars.janusPatchF[firstIPC][i] -= modulus;
                loopVars.janusPatchF[secndIPC][i] += modulus;
            }
        }
    }
}
