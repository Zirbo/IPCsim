#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iomanip>
#include "IPCsimulation.hpp"



void IPCsimulation::printPotentialsToFile(int potentialPrintingStep) {
    if(isJanusSimulation) {
        printPotentialsToFileJanus(potentialPrintingStep);
        return;
    }

    // interactionRange and forceAndEnergySamplingStep are both scaled by simulationBoxSide, so their ratio is right
    const size_t potentialRangeSamplingSize = size_t( interactionRange/forceAndEnergySamplingStep ) + 1;
    const size_t numberOfPrints = potentialRangeSamplingSize/potentialPrintingStep;

    const std::string dirName("potentials_for_lammps/");
    const std::string extension(".table");
    std::string interactionType [6];
    interactionType[0] = "BB";      interactionType[1] = "Bs1";     interactionType[2] = "Bs2";
    interactionType[3] = "s1s2";    interactionType[4] = "s1s1";    interactionType[5] = "s2s2";

    for (int type = 0; type < 6; ++type) {
        std::string fileName = dirName + interactionType[type] + extension;
        std::ofstream potentialOutputFile(fileName.c_str());
        potentialOutputFile << "# potentials for lammps\n\n" << interactionType[type] << "\nN " << numberOfPrints << "\n\n";
        potentialOutputFile << std::scientific << std::setprecision(6);

        int printCounter = 0;
        for ( size_t i = potentialPrintingStep; i < potentialRangeSamplingSize; i += potentialPrintingStep) {
            const double r = i*forceAndEnergySamplingStep*simulationBoxSide;
            printCounter++;
            potentialOutputFile << printCounter << "\t" << r << "\t";
            if( type == 0) {
                potentialOutputFile << uBB[i]*r << "\t" << fBB[i]*r << "\n";
            } else if ( type == 1) {
                potentialOutputFile << uBs1[i]*r << "\t" << fBs1[i]*r << "\n";
            } else if ( type == 2) {
                potentialOutputFile << uBs2[i]*r << "\t" << fBs2[i]*r << "\n";
            } else if ( type == 3) {
                potentialOutputFile << us1s2[i]*r << "\t" << fs1s2[i]*r << "\n";
            } else if ( type == 4) {
                potentialOutputFile << us1s1[i]*r << "\t" << fs1s1[i]*r << "\n";
            } else if ( type == 5) {
                potentialOutputFile << us2s2[i]*r << "\t" << fs2s2[i]*r << "\n";
            }
        }
        potentialOutputFile.close();
    }
}


void IPCsimulation::computeVerletHalfStepForIPC(IPC & ipc) {
    double dxNew[3], dxOld[3];
    for (int i: {0, 1, 2}) {
        // compute and store the old separation, so that we no longer need the old positions.
        dxOld[i] = ipc.firstPatch.x[i] - ipc.secndPatch.x[i];
        relativePBC(dxOld[i]);

        // compute the half step velocities from the effective forces of the last step
        ipc.firstPatch.v[i] += ipc.eFp1[i]*halfDtFirstPatchInverseMass;
        ipc.secndPatch.v[i] += ipc.eFp2[i]*halfDtSecndPatchInverseMass;

        // compute the new positions from the half step velocities
        ipc.firstPatch.x[i] += ipc.firstPatch.v[i]*dt;
        absolutePBC(ipc.firstPatch.x[i]);
        ipc.secndPatch.x[i] += ipc.secndPatch.v[i]*dt;
        absolutePBC(ipc.secndPatch.x[i]);

        // compute the new separation between the two patches
        dxNew[i] = ipc.firstPatch.x[i] - ipc.secndPatch.x[i] ;
        relativePBC(dxNew[i]);
    }
    // compute the (squared) violation of the constraint
    double diff = (dxNew[0]*dxNew[0] + dxNew[1]*dxNew[1] + dxNew[2]*dxNew[2]) - squaredPatchDistance;

    // correct the positions and the velocities until the violation is less than the tollerance
    while( std::fabs(diff) > tollerance*squaredPatchDistance )
    {
        double g = diff/( 2.*(dxOld[0]*dxNew[0] + dxOld[1]*dxNew[1] + dxOld[2]*dxNew[2]) * alpha_sum*dt );

        for (int i: {0, 1, 2}) {
            double DXi = g*dxOld[i];

            ipc.firstPatch.v[i] -= alpha_1*DXi;
            ipc.secndPatch.v[i] += alpha_2*DXi;

            DXi *= dt;

            ipc.firstPatch.x[i] -= DXi*alpha_1;
            absolutePBC(ipc.firstPatch.x[i]);
            ipc.secndPatch.x[i]  += DXi*alpha_2;
            absolutePBC(ipc.secndPatch.x[i] );

            dxNew[i] = ipc.firstPatch.x[i] - ipc.secndPatch.x[i];
            relativePBC(dxNew[i]);
        }

        diff = (dxNew[0]*dxNew[0] + dxNew[1]*dxNew[1] + dxNew[2]*dxNew[2]) - squaredPatchDistance;
    }

    for (int i: {0, 1, 2}) {
        ipc.ipcCenter.x[i] = ipc.secndPatch.x[i] + dxNew[i]*secndPatchEccentricity*inversePatchDistance;
        absolutePBC(ipc.ipcCenter.x[i]);
    }
}


void IPCsimulation::finishVerletStepForIPC(IPC & ipc) {
    double dx[3], dv[3];
    for (int i: {0, 1, 2}) {
        if (std::fabs(ipc.eFp1[i]) > 1e5) {
            return;
        }
        // compute the the final velocities from the new effective forces
        ipc.firstPatch.v[i] += ipc.eFp1[i]*halfDtFirstPatchInverseMass;
        ipc.secndPatch.v[i] += ipc.eFp2[i]*halfDtSecndPatchInverseMass;

        // compute the patch-patch distance
        dx[i] = ipc.firstPatch.x[i] - ipc.secndPatch.x[i];
        relativePBC(dx[i]);
        dv[i] = ipc.firstPatch.v[i] - ipc.secndPatch.v[i];
    }
    // check how much the constraints are being violated
    double k = (dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2])*inverseAlpha_sumSquaredPatchDistance;
    while( std::fabs(k) > tollerance ) {
        // compute and apply corrections
        for (int i: {0, 1, 2}) {
            const double DVi = k*dx[i];
            ipc.firstPatch.v[i] -= DVi*alpha_1;
            ipc.secndPatch.v[i] += DVi*alpha_2;
            dv[i] = ipc.firstPatch.v[i] - ipc.secndPatch.v[i];
        }
        // recompute the violation of the constraints
        k = (dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2])*inverseAlpha_sumSquaredPatchDistance;
    }

    for (int i: {0, 1, 2}) {
        ipc.ipcCenter.v[i] = (ipc.firstPatch.v[i]*secndPatchEccentricity + ipc.secndPatch.v[i]*firstPatchEccentricity)*inversePatchDistance;
    }
}



void IPCsimulation::initializeNewConfiguration(int N1) {
    if(isJanusSimulation) {
        initializeNewJanusConfiguration(N1);
        return;
    }

    particles.resize(nIPCs);
    RandomNumberGenerator rand;

    int N2 = N1*N1;
    int N3 = N2*N1;
    double vel_scaling = std::sqrt(1.6*initialTemperature)/simulationBoxSide;

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
        generateRandomOrientation(ipcAxis,rand);
        generateRandomOrientation(ipcOrthogonalAxis, rand);
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
    if(isJanusSimulation) {
        restorePreviousJanusConfiguration();
        return;
    }

    char unusedPatchName;
    double unusedTime;
    std::ifstream startingConfigurationFile("startingstate.xyz");
    if(startingConfigurationFile.fail()) {
        std::cerr << "File startingstate.xyz could not be opened. Aborting.\n";
        exit(1);
    }
    startingConfigurationFile >> nIPCs >> simulationBoxSide >> unusedTime;
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

    if(isJanusSimulation) {
        computeFreeJanusForces();
        return;
    }

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
