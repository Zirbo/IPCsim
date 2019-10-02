#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iomanip>
#include "IPCsimulation.hpp"



void IPCsimulation::printPotentialsToFileJanus(int potentialPrintingStep) {

    // interactionRange and forceAndEnergySamplingStep are both scaled by simulationBoxSide, so their ratio is right
    const size_t potentialRangeSamplingSize = size_t( interactionRange/forceAndEnergySamplingStep ) + 1;
    const size_t numberOfPrints = potentialRangeSamplingSize/potentialPrintingStep;

    const std::string dirName("potentials_for_lammps/");
    const std::string extension(".table");
    std::string interactionType [3];
    interactionType[0] = "BB";      interactionType[1] = "Bs";      interactionType[2] = "ss";

    for (int type = 0; type < 3; ++type) {
        std::string fileName = dirName + interactionType[type] + extension;
        std::ofstream potentialOutputFile(fileName);
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
                potentialOutputFile << us1s2[i]*r << "\t" << fs1s2[i]*r << "\n";
            }
        }
        potentialOutputFile.close();
    }
}


/*****************************************************************************************/
void IPCsimulation::computeVerletHalfStepForJanusIPC(JanusIPC & ipc) {
    double xc[3], xp[3];
    double dxNew[3];
    for (int i: {0, 1, 2}) {
        // compute the half step velocities from the effective forces of the last step
        ipc.ipcCenter.v[i] += ipc.ipcCenter.F[i]*(.5*dt*ipcCenterInverseMass);
        ipc.janusPatch.v[i] += ipc.janusPatch.F[i]*(.5*dt*firstPatchInverseMass);

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
    double diff = (dxNew[0]*dxNew[0] + dxNew[1]*dxNew[1] + dxNew[2]*dxNew[2]) - squaredPatchDistance;

    // correct the positions and the velocities until the violation is less than the tollerance
    while( std::fabs(diff) > tollerance*squaredPatchDistance )
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

        diff = (dxNew[0]*dxNew[0] + dxNew[1]*dxNew[1] + dxNew[2]*dxNew[2]) - squaredPatchDistance;
    }

    for (int i: {0, 1, 2}) {
        ipc.ipcCenter.x[i] = xc[i];
        ipc.janusPatch.x[i] = xp[i];
    }
}


void IPCsimulation::finishVerletStepForJanusIPC(JanusIPC & ipc) {
    double vc[3], vp[3], dx[3], dv[3];
    for (int i: {0, 1, 2}) {
        // compute the the final velocities from the new effective forces
        vc[i] = ipc.ipcCenter.v[i] + ipc.ipcCenter.F[i]*(.5*dt*ipcCenterInverseMass);
        vp[i] = ipc.janusPatch.v[i] + ipc.janusPatch.F[i]*(.5*dt*firstPatchInverseMass);

        // compute the patch-patch distance
        dx[i] = ipc.ipcCenter.x[i] - ipc.janusPatch.x[i];
        relativePBC(dx[i]);
        dv[i] = vc[i] - vp[i];
    }
    // check how much the constraints are being violated
    double k = (dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2]) / squaredPatchDistance;
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
        k = (dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2]) / squaredPatchDistance;
    }

    for (int i: {0, 1, 2}) {
        ipc.ipcCenter.v[i] = vc[i];
        ipc.janusPatch.v[i] = vp[i];
    }
}




void IPCsimulation::initializeNewJanusConfiguration(int N1) {
    janusParticles.resize(nIPCs);
    RandomNumberGenerator rand;

    int N2 = N1*N1;
    int N3 = N2*N1;

    // scaling: sqrt(2kT/mPI) comes from boltzmann average of |v_x|
    //double vel_scaling = std::sqrt(2.*desiredTemperature/3.1415)/simulationBoxSide;

    double vel_scaling = std::sqrt(1.6*initialTemperature)/simulationBoxSide;
    // initialize IPC positions
    for(int i=0;i<N3;i++)
    {
      // FCC is obtained as 4 intersecating SC
        janusParticles[i].number = i;
        janusParticles[i].type = 'C';
        janusParticles[i].ipcCenter.x[0] = (i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(janusParticles[i].ipcCenter.x[0]);
        janusParticles[i].ipcCenter.x[1] = ((i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(janusParticles[i].ipcCenter.x[1]);
        janusParticles[i].ipcCenter.x[2] = (i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(janusParticles[i].ipcCenter.x[2]);

        janusParticles[i+N3].number = i+N3;
        janusParticles[i+N3].type = 'C';
        janusParticles[i+N3].ipcCenter.x[0] = (.5 + i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(janusParticles[i+N3].ipcCenter.x[0]);
        janusParticles[i+N3].ipcCenter.x[1] = (.5 + (i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(janusParticles[i+N3].ipcCenter.x[1]);
        janusParticles[i+N3].ipcCenter.x[2] = (i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(janusParticles[i+N3].ipcCenter.x[2]);

        janusParticles[i+N3+N3].number = i+N3+N3;
        janusParticles[i+N3+N3].type = 'C';
        janusParticles[i+N3+N3].ipcCenter.x[0] = (i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(janusParticles[i+N3+N3].ipcCenter.x[0]);
        janusParticles[i+N3+N3].ipcCenter.x[1] = (.5 + (i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(janusParticles[i+N3+N3].ipcCenter.x[1]);
        janusParticles[i+N3+N3].ipcCenter.x[2] = (.5 + i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(janusParticles[i+N3+N3].ipcCenter.x[2]);

        janusParticles[i+N3+N3+N3].number = i+N3+N3+N3;
        janusParticles[i+N3+N3+N3].type = 'C';
        janusParticles[i+N3+N3+N3].ipcCenter.x[0] = (.5 + i%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(janusParticles[i+N3+N3+N3].ipcCenter.x[0]);
        janusParticles[i+N3+N3+N3].ipcCenter.x[1] = ((i/N1)%N1 + .1*rand.getRandom55())/N1;
        absolutePBC(janusParticles[i+N3+N3+N3].ipcCenter.x[1]);
        janusParticles[i+N3+N3+N3].ipcCenter.x[2] = (.5 + i/N2 + .1*rand.getRandom55()) /N1;
        absolutePBC(janusParticles[i+N3+N3+N3].ipcCenter.x[2]);

        // starting from random but ONLY TRANSLATIONAL speeds, for compatibility with rattle
        for (int j: {0, 1, 2}) {
           janusParticles[i].ipcCenter.v[j]          = rand.getRandom11()*vel_scaling;
           janusParticles[i+N3].ipcCenter.v[j]       = rand.getRandom11()*vel_scaling;
           janusParticles[i+N3+N3].ipcCenter.v[j]    = rand.getRandom11()*vel_scaling;
           janusParticles[i+N3+N3+N3].ipcCenter.v[j] = rand.getRandom11()*vel_scaling;
        }
    }
    // initialize patches positions
    for(JanusIPC &ipc: janusParticles) {
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

            ipc.janusPatch.x[i] = ipc.ipcCenter.x[i] + ipcAxis[i]*patchDistance;
            absolutePBC(ipc.janusPatch.x[i]);

            double temp = rand.getRandom11()*ipcOrthogonalAxis[i]*vel_scaling;
            ipc.janusPatch.v[i] = ipc.ipcCenter.v[i] + temp;
        }
    }
}

void IPCsimulation::restorePreviousJanusConfiguration() {
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

    for (JanusIPC &ipc: janusParticles) {
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



void IPCsimulation::computeFreeJanusForces() {
    // reset all forces
    for(JanusIPC &ipc: janusParticles) {
        for (int i: {0, 1, 2}) {
            ipc.ipcCenter.F[i] = 0.;
            ipc.janusPatch.F[i] = 0.;
        }
    }
    potentialEnergy = 0.0;  squaredMinimumDistanceBetweenParticles = 1.;

    #pragma omp parallel
    {
        loopVariablesJanus loopVars(nIPCs);

        #pragma omp for
        for(int m=0; m<cells.getNumberofCells(); ++m)  // loop over all cells
        {
            const std::list<int> & ipcsInCell = cells.getIPCsInCell(m);
            const std::list<int> ipcsInNeighbouringCells = cells.getIPCsInNeighbouringCells(m);
            for(auto ipc = ipcsInCell.cbegin(); ipc != ipcsInCell.cend(); ++ipc) {
                computeInteractionsWithJanusIPCsInTheSameCell(ipc, ipcsInCell, loopVars);
                computeInteractionsWithJanusIPCsInNeighbouringCells(ipc, ipcsInNeighbouringCells, loopVars);
            }
        }
        #pragma omp critical
        {
            for (JanusIPC &ipc: janusParticles) {
                for (int i: {0, 1, 2}) {
                    const size_t j = ipc.number;
                    ipc.ipcCenter.F[i]  += loopVars.ipcCenterF[j][i];
                    ipc.janusPatch.F[i] += loopVars.janusPatchF[j][i];
                }
            }
            potentialEnergy += loopVars.U;
            squaredMinimumDistanceBetweenParticles += loopVars.minimumSquaredDistance;
        }

        for(JanusIPC &ipc: janusParticles) {
            for (int i: {0, 1, 2}) {
                if(isFieldEnabled) {
                    ipc.ipcCenter.F[i] += externalFieldIpcCenter[i];
                    ipc.janusPatch.F[i] += externalFieldFirstPatch[i];
                }
            }
        }
    }
}

void IPCsimulation::computeInteractionsWithJanusIPCsInNeighbouringCells(std::list<int>::const_iterator loc, std::list<int> const& ipcsInNeighbouringCells, loopVariablesJanus & loopVars) {
    for( auto ext = ipcsInNeighbouringCells.cbegin(); ext != ipcsInNeighbouringCells.cend(); ++ext) {
        computeInteractionsBetweenTwoJanusIPCs(*loc, *ext, loopVars);
    }
}



void IPCsimulation::computeInteractionsWithJanusIPCsInTheSameCell(std::list<int>::const_iterator loc, std::list<int> const& ipcsInCurrentCell, loopVariablesJanus &loopVars) {
    // starts from loc+1 which is like summing over i > j inside the cell
    for(std::list<int>::const_iterator ins = std::next(loc); ins != ipcsInCurrentCell.cend(); ++ins) {
        computeInteractionsBetweenTwoJanusIPCs(*loc, *ins, loopVars);
    }
}

void IPCsimulation::computeInteractionsBetweenTwoJanusIPCs(const int firstIPC, const int secndIPC, loopVariablesJanus &loopVars) {

    JanusIPC const& first = janusParticles[firstIPC];
    JanusIPC const& secnd = janusParticles[secndIPC];

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
            loopVars.U += uBs1[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fBs1[dist]*siteSiteSeparation[0][i];
                loopVars.ipcCenterF[firstIPC][i] -= modulus;
                loopVars.janusPatchF[secndIPC][i] += modulus;
            }
        } else if (j == 1) { // patch - center
            loopVars.U += uBs1[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fBs1[dist]*siteSiteSeparation[1][i];
                loopVars.janusPatchF[firstIPC][i] -= modulus;
                loopVars.ipcCenterF[secndIPC][i] += modulus;
            }
        } else if (j == 2) { // patch - patch
            loopVars.U += us1s1[dist];
            for (int i: {0, 1, 2}) {
                const double modulus = fs1s1[dist]*siteSiteSeparation[2][i];
                loopVars.janusPatchF[firstIPC][i] -= modulus;
                loopVars.janusPatchF[secndIPC][i] += modulus;
            }
        }
    }
}
