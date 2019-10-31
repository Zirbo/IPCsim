#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "IPCsimulation.hpp"


void IPCsimulation::printRawSiteSitePotentials(const int potentialPrintingStep) {
    std::ofstream potentialOutputFile("potentials/site-site-potentials.txt");
    potentialOutputFile << "# r\t\tHS\t\tBB\t\tBs1\t\tBs2\t\ts1s1\t\ts1s2\t\ts2s2\n";
    potentialOutputFile << std::scientific << std::setprecision(6);

    std::ofstream forcesOutputFile("potentials/site-site-forces.txt");
    forcesOutputFile << "# r\t\tHS\t\tBB\t\tBs1\t\tBs2\t\ts1s1\t\ts1s2\t\ts2s2\n";
    forcesOutputFile << std::scientific << std::setprecision(6);

    // interactionRange and forceAndEnergySamplingStep are both scaled by simulationBoxSide, so their ratio is right
    const size_t potentialRangeSamplingSize = size_t( interactionRange/forceAndEnergySamplingStep ) + 1;

    for ( size_t i = potentialPrintingStep; i < potentialRangeSamplingSize; i += potentialPrintingStep) {
        const double r = i*forceAndEnergySamplingStep*simulationBoxSide;
        potentialOutputFile << r << "\t" << uHS[i] << "\t"
                            << uBB[i] << "\t" << uBs1[i] << "\t" << uBs2[i] << "\t"
                            << us1s2[i] << "\t" << us1s1[i] << "\t" << us2s2[i] << "\n";
        forcesOutputFile << r << "\t" << fHS[i] << "\t"
                         << fBB[i]*r << "\t" << fBs1[i]*r << "\t" << fBs2[i]*r << "\t"
                         << fs1s2[i]*r << "\t" << fs1s1[i]*r << "\t" << fs2s2[i]*r << "\n";
    }
    potentialOutputFile.close();
}


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
                    potentialOutputFile << uHS[i] + uBB[i] << "\t" << (fHS[i] + fBB[i])*r << "\n";
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
                    printPotential = ( uHS[i] + uBB[i] > cutoffValue )? cutoffValue : uHS[i] + uBB[i];
                    printForce = ((fHS[i] + fBB[i])*r < -cutoffValue )? -cutoffValue : (fHS[i] + fBB[i])*r;
                } else if ( type == 1) {
                    printPotential = ( uBs1[i] > cutoffValue )? cutoffValue : uBs1[i];
                    printForce = (fBs1[i]*r < -cutoffValue )? -cutoffValue : fBs1[i]*r;
                } else if ( type == 2) {
                    printPotential = ( uBs2[i] > cutoffValue )? cutoffValue : uBs2[i];
                    printForce = (fBs2[i]*r < -cutoffValue )? -cutoffValue : fBs2[i]*r;
                } else if ( type == 3) {
                    printPotential = ( us1s2[i] > cutoffValue )? cutoffValue : us1s2[i];
                    printForce = (fs1s2[i]*r < -cutoffValue )? -cutoffValue : fs1s2[i]*r;
                } else if ( type == 4) {
                    printPotential = ( us1s1[i] > cutoffValue )? cutoffValue : us1s1[i];
                    printForce = (fs1s1[i]*r < -cutoffValue )? -cutoffValue : fs1s1[i]*r;
                } else if ( type == 5) {
                    printPotential = ( us2s2[i] > cutoffValue )? cutoffValue : us2s2[i];
                    printForce = (fs2s2[i]*r < -cutoffValue )? -cutoffValue : fs2s2[i]*r;
                }
                potentialOutputFile << printPotential << "\t" << printForce << "\n";
            }
        }
        potentialOutputFile.close();
    }

}


void IPCsimulation::printPotentialsToFileForVisualization(const int potentialPrintingStep) {
    particles.resize(2);
    // compute the four possible orientation with two asymmetric patches
    IPC upwards, downwards, sidewaysLeft, sidewaysRight;

    upwards.type = 'u';
    upwards.ipcCenter  = {0.5, 0.5, 0.5};
    upwards.firstPatch = {0.5, 0.5, 0.5 + firstPatchEccentricity};
    upwards.secndPatch = {0.5, 0.5, 0.5 - secndPatchEccentricity};

    downwards.type = 'd';
    downwards.ipcCenter  = {0.5, 0.5, 0.5};
    downwards.firstPatch = {0.5, 0.5, 0.5 - firstPatchEccentricity};
    downwards.secndPatch = {0.5, 0.5, 0.5 + secndPatchEccentricity};

    sidewaysLeft.type = 'l';
    sidewaysLeft.ipcCenter  = {0.5, 0.5, 0.5};
    sidewaysLeft.firstPatch = {0.5 + firstPatchEccentricity, 0.5, 0.5};
    sidewaysLeft.secndPatch = {0.5 - secndPatchEccentricity, 0.5, 0.5};

    sidewaysRight.type = 'r';
    sidewaysRight.ipcCenter  = {0.5, 0.5, 0.5};
    sidewaysRight.firstPatch = {0.5 - firstPatchEccentricity, 0.5, 0.5};
    sidewaysRight.secndPatch = {0.5 + secndPatchEccentricity, 0.5, 0.5};

    // the couples of independent orientations are: uu, ud, ul, ur, ll, lr, rr (because dd=-uu, dl=-ul, dr=-ur)
    typedef std::pair<IPC,IPC> IPCcouple;
    for( IPCcouple const& couple: { IPCcouple(upwards, upwards), IPCcouple(upwards, downwards),
                        IPCcouple(upwards, sidewaysLeft), IPCcouple(upwards, sidewaysRight),
                        IPCcouple(sidewaysLeft, sidewaysLeft), IPCcouple(sidewaysLeft, sidewaysRight),
                        IPCcouple(sidewaysRight, sidewaysLeft) }) {
        particles[0] = couple.first;
        particles[1] = couple.second;
        printPotentialsToFileForVisualizationSingleOrientation(potentialPrintingStep);
    }
}

void IPCsimulation::printPotentialsToFileForVisualizationSingleOrientation(const int potentialPrintingStep) {
    // open two output files
    const std::string dirName("potentials/");
    const std::string extension(".dat");
    std::string fileName = dirName + "emanuela_r_" + particles[0].type + particles[1].type + extension;
    std::ofstream radialOutputFile(fileName);
    fileName = dirName + "countour_plot_" + particles[0].type + particles[1].type + extension;
    std::ofstream contourPlotOutputFile(fileName);

    radialOutputFile << "# potentials plotted as f(r) in Emanuela's six independent orientations\n";
    radialOutputFile << "# r          U(r)           U(r)*e_min (i.e. not renormalized)\n";
    radialOutputFile << std::scientific << std::setprecision(6);

    contourPlotOutputFile << "# potentials plotted as f(x,z) on a half-plane in Emanuela's six indepentent orientations\n";
    contourPlotOutputFile << "# x          z            U(x,z)           U(x,z)*e_min (i.e. not renormalized)\n";
    contourPlotOutputFile << std::scientific << std::setprecision(6);


    const double dr = forceAndEnergySamplingStep*potentialPrintingStep;
    const double atContact = 1./simulationBoxSide;

    const IPC initialStateProbeIpc = particles[1];

    particles[1].ipcCenter.x[2]  -= interactionRange;
    particles[1].firstPatch.x[2] -= interactionRange;
    particles[1].secndPatch.x[2] -= interactionRange;

    for (double z = -interactionRange; z <= interactionRange; z += dr) {
        for (double x = 0; x <= interactionRange; x += dr) {
            loopVariables loopVars(2);
            computeInteractionsBetweenTwoIPCs(0, 1, loopVars);

            double contourPotential = (loopVars.U < 5.)? loopVars.U : 5.;
            if (contourPotential < -1.) {
                std::cerr << std::scientific << std::setprecision(6);
                std::cerr << contourPotential << " < -1.0!!!!\n";
                contourPotential = -1.;
            }
            contourPlotOutputFile << x*simulationBoxSide << "\t" << z*simulationBoxSide << "\t" << contourPotential << "\t" << contourPotential*e_min << "\n";

            // if z = 0 and 1. < x < interaction range, print the six site-site potentials
            if (std::fabs(z) < dr*tollerance) {
                if (x > atContact - dr) {
                    radialOutputFile << x*simulationBoxSide << "\t" << loopVars.U << "\t" << loopVars.U*e_min << "\n";
                }
            }

            // update x-positions
            particles[1].ipcCenter.x[0]  += dr;
            particles[1].firstPatch.x[0] += dr;
            particles[1].secndPatch.x[0] += dr;
        }

        // reset x-positions, update z-positions
        particles[1].ipcCenter.x[0]  = initialStateProbeIpc.ipcCenter.x[0];
        particles[1].firstPatch.x[0] = initialStateProbeIpc.firstPatch.x[0];
        particles[1].secndPatch.x[0] = initialStateProbeIpc.secndPatch.x[0];

        particles[1].ipcCenter.x[2]  += dr;
        particles[1].firstPatch.x[2] += dr;
        particles[1].secndPatch.x[2] += dr;
    }
}
