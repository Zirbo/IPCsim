#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "IPCsimulation.hpp"

void IPCsimulation::printPotentialsToFileLAMMPS(const int potentialPrintingStep, const int cutoffValue) {
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
        std::ofstream potentialOutputFile(fileName);
        potentialOutputFile << "# potentials for lammps\n\n" << interactionType[type] << "\nN " << numberOfPrints << "\n\n";
        potentialOutputFile << std::scientific << std::setprecision(6);

        int printCounter = 0;
        for ( size_t i = potentialPrintingStep; i < potentialRangeSamplingSize; i += potentialPrintingStep) {
            const double r = i*forceAndEnergySamplingStep*simulationBoxSide;
            printCounter++;
            potentialOutputFile << printCounter << "\t" << r << "\t";

            if (cutoffValue < 0) {
                if( type == 0) {
                    potentialOutputFile << uBB[i] << "\t" << fBB[i]*r << "\n";
                } else if ( type == 1) {
                    potentialOutputFile << uBs1[i] << "\t" << fBs1[i]*r << "\n";
                } else if ( type == 2) {
                    potentialOutputFile << uBs2[i] << "\t" << fBs2[i]*r << "\n";
                } else if ( type == 3) {
                    potentialOutputFile << us1s2[i] << "\t" << fs1s2[i]*r << "\n";
                } else if ( type == 4) {
                    potentialOutputFile << us1s1[i] << "\t" << fs1s1[i]*r << "\n";
                } else if ( type == 5) {
                    potentialOutputFile << us2s2[i] << "\t" << fs2s2[i]*r << "\n";
                }
            } else {
                double printPotential, printForce;
                if( type == 0) {
                    printPotential = ( uBB[i] > cutoffValue )? cutoffValue : uBB[i];
                    printForce = (fBB[i]*r < -cutoffValue )? -cutoffValue : fBB[i];
                } else if ( type == 1) {
                    printPotential = ( uBs1[i] > cutoffValue )? cutoffValue : uBs1[i];
                    printForce = (fBs1[i]*r < -cutoffValue )? -cutoffValue : fBs1[i];
                } else if ( type == 2) {
                    printPotential = ( uBs2[i] > cutoffValue )? cutoffValue : uBs2[i];
                    printForce = (fBs2[i]*r < -cutoffValue )? -cutoffValue : fBs2[i];
                } else if ( type == 3) {
                    printPotential = ( us1s2[i] > cutoffValue )? cutoffValue : us1s2[i];
                    printForce = (fs1s2[i]*r < -cutoffValue )? -cutoffValue : fs1s2[i];
                } else if ( type == 4) {
                    printPotential = ( us1s1[i] > cutoffValue )? cutoffValue : us1s1[i];
                    printForce = (fs1s1[i]*r < -cutoffValue )? -cutoffValue : fs1s1[i];
                } else if ( type == 5) {
                    printPotential = ( us2s2[i] > cutoffValue )? cutoffValue : us2s2[i];
                    printForce = (fs2s2[i]*r < -cutoffValue )? -cutoffValue : fs2s2[i];
                }
                potentialOutputFile << printPotential << "\t" << printForce << "\n";
            }
        }
        potentialOutputFile.close();
    }

}


void IPCsimulation::printPotentialsToFileEmanuela6r(const int potentialPrintingStep) {
    particles.resize(2);
    // compute the four possible orientation with two asymmetric patches
    IPC upwards, downwards, sidewaysLeft, sidewaysRight;

    upwards.type = 'u';
    upwards.ipcCenter  = {0., 0., 0.};
    upwards.firstPatch = {0., 0., firstPatchEccentricity};
    upwards.secndPatch = {0., 0., -secndPatchEccentricity};

    downwards.type = 'd';
    downwards.ipcCenter  = {0., 0., 0.};
    downwards.firstPatch = {0., 0., -firstPatchEccentricity};
    downwards.secndPatch = {0., 0., secndPatchEccentricity};

    sidewaysLeft.type = 'l';
    sidewaysLeft.ipcCenter  = {0., 0., 0.};
    sidewaysLeft.firstPatch = {firstPatchEccentricity, 0., 0.};
    sidewaysLeft.secndPatch = {-secndPatchEccentricity, 0., 0.};

    sidewaysRight.type = 'r';
    sidewaysRight.ipcCenter  = {0., 0., 0.};
    sidewaysRight.firstPatch = {-firstPatchEccentricity, 0., 0};
    sidewaysRight.secndPatch = {secndPatchEccentricity, 0., 0.};

    // the couples of orientations that are not equal by symmetry are: uu-00 ud-01 ul-02 ur-03 ll-22 lr-23 (because uu=-dd, dl=-ul, dr=-ur, rr=-ll)
    typedef std::pair<IPC,IPC> IPCcouple;
    for( IPCcouple const& couple: { IPCcouple(upwards, upwards), IPCcouple(upwards, downwards),
                        IPCcouple(upwards, sidewaysLeft), IPCcouple(upwards, sidewaysRight),
                        IPCcouple(sidewaysLeft, sidewaysLeft), IPCcouple(sidewaysLeft, sidewaysRight) }) {
        particles[0] = couple.first;
        particles[1] = couple.second;
        printPotentialsToFileEmanuela6rSingleOrientation(potentialPrintingStep);
    }
}

void IPCsimulation::printPotentialsToFileEmanuela6rSingleOrientation(const int potentialPrintingStep) {

    const std::string dirName("potentials_emanuela6r/");
    const std::string extension(".table");
    const std::string fileName = dirName + particles[0].type + particles[1].type + extension;
    //std::stringstream fileName;
    //fileName << "potentials_emanuela6r/" << std::string(1, particles[0].type) << std::string(1, particles[1].type) << ".table";
    //const std::string fileName = "potentials_emanuela6r/" + particles[0].type + particles[1].type + ".table";
    std::ofstream potentialOutputFile(fileName);
    potentialOutputFile << "# potentials plotted as f(r) in Emanuela's six indipentent orientations\n";
    potentialOutputFile << "# r          U(r)           U(r)*e_min (not renormalized)\n";
    potentialOutputFile << std::scientific << std::setprecision(6);

    const double dr = forceAndEnergySamplingStep*potentialPrintingStep;
    const double atContact = 1./simulationBoxSide;
    particles[1].ipcCenter.x[0]  += atContact;
    particles[1].firstPatch.x[0] += atContact;
    particles[1].secndPatch.x[0] += atContact;

    for (double r = atContact; r <= interactionRange; r += dr) {
        loopVariables loopVars(2);
        computeInteractionsBetweenTwoIPCs(0, 1, loopVars);

        potentialOutputFile << r*simulationBoxSide << "\t" << loopVars.U << "\t" << loopVars.U*e_min << "\n";

        particles[1].ipcCenter.x[0]  += dr;
        particles[1].firstPatch.x[0] += dr;
        particles[1].secndPatch.x[0] += dr;
    }
}



void IPCsimulation::printPotentialsToFileContourPlot(const int potentialPrintingStep) {
    // not yet implemented
}
