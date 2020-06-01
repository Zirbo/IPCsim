#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

void error(std::string const & errorMessage) {
    std::cerr << errorMessage;
    exit(1);
}

class PotentialForLammps {
public:
    PotentialForLammps(std::string const& inputFileName);
    void computePotential();
    void printPotentialsToFile(std::string const& outputDirName);

private:
    double e_BB, e_Bs1, e_Bs2, e_s1s1, e_s1s2, e_s2s2, e_min;
    double firstPatchEccentricity, firstPatchRadius;
    double secndPatchEccentricity, secndPatchRadius;
    double ipcRadius, interactionRange;
    double fakeHSdiameter, fakeHScoefficient, fakeHSexponent;

    double samplingStep, cutoffValue;

    std::vector<double> uHS, uBB, uBs1, uBs2, us1s2, us1s1, us2s2;
    std::vector<double> fHS, fBB, fBs1, fBs2, fs1s2, fs1s1, fs2s2;

    double computeOmega(double Ra, double Rb, double rab);
    double computeOmegaRadialDerivative(double Ra, double Rb, double rab);
};

PotentialForLammps::PotentialForLammps(const std::string &inputFileName) {
    // read input.in file
    std::ifstream inputFile(inputFileName);
    if(inputFile.fail())
        error("File input.in could not be opened. Aborting.\n");
    inputFile >> e_BB >> e_Bs1 >> e_Bs2;
    inputFile >> e_s1s1 >> e_s2s2 >> e_s1s2;
    inputFile >> e_min;
    inputFile >> firstPatchEccentricity >> firstPatchRadius;
    inputFile >> secndPatchEccentricity >> secndPatchRadius;
    inputFile >> fakeHScoefficient >> fakeHSexponent;
    inputFile >> samplingStep >> cutoffValue;
    inputFile.close();

    // patch geometry integrity check
    if ( std::abs( ( firstPatchEccentricity + firstPatchRadius ) - ( secndPatchEccentricity + secndPatchRadius ) ) >= 1e-10 )
        error("eccentricities and radii are not consistent!\n");

    fakeHSdiameter = 1.0;
    ipcRadius = firstPatchEccentricity + firstPatchRadius;
    interactionRange = 2*ipcRadius;
}

double PotentialForLammps::computeOmega(double Ra, double Rb, double rab) {
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

double PotentialForLammps::computeOmegaRadialDerivative(double Ra, double Rb, double rab) {
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

void PotentialForLammps::computePotential()
{
    const size_t potentialSteps = size_t( interactionRange/samplingStep ) + 2;

    uHS.resize(potentialSteps);
    uBB.resize(potentialSteps);
    uBs1.resize(potentialSteps);
    uBs2.resize(potentialSteps);
    us1s2.resize(potentialSteps);
    us1s1.resize(potentialSteps);
    us2s2.resize(potentialSteps);

    fHS.resize(potentialSteps);
    fBB.resize(potentialSteps);
    fBs1.resize(potentialSteps);
    fBs2.resize(potentialSteps);
    fs1s2.resize(potentialSteps);
    fs1s1.resize(potentialSteps);
    fs2s2.resize(potentialSteps);

    for ( size_t i = 0; i < potentialSteps; ++i)
    {
        const double r = i*samplingStep;
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

        if ( r <= fakeHSdiameter )
        {
            // setting up a Fake Hard Sphere Core
            double rm = pow(r, -fakeHSexponent);
            uHS[i]   += fakeHScoefficient*((rm-2.)*rm+1.);
            fHS[i]   += -2.*fakeHSexponent*fakeHScoefficient*(rm-1.)*rm;
        }
    }
}

void PotentialForLammps::printPotentialsToFile(const std::string &outputDirName) {
    // create output directory
    const std::string makedir = "mkdir -p " + outputDirName;
    if(system(makedir.c_str()) != 0)
        error("Problem while creating the directory.\n");

    // prepare strings for defining file names
    const std::string extension(".table");
    std::string interactionType [6];
    interactionType[0] = "/BB";      interactionType[1] = "/Bs1";     interactionType[2] = "/Bs2";
    interactionType[3] = "/s1s2";    interactionType[4] = "/s1s1";    interactionType[5] = "/s2s2";

    const size_t potentialSteps = uHS.size();

    for (int type = 0; type < 6; ++type) {
        // create the output file
        std::string fileName = outputDirName + interactionType[type] + extension;
        std::ofstream potentialOutputFile(fileName);
        potentialOutputFile << "# potentials for lammps\n\n" << interactionType[type] << "\nN " << potentialSteps << "\n\n";
        potentialOutputFile << std::scientific << std::setprecision(6);

        for ( size_t i = 1; i < potentialSteps; ++i) {
            const double r = i*samplingStep;
            potentialOutputFile << i << "\t" << r << "\t";

            // compute force and potential depending on type and cutoff
            double printPotential, printForce;
            if( type == 0) {
                printPotential = ( uHS[i] + uBB[i] > cutoffValue )? cutoffValue : uHS[i] + uBB[i];
                printForce = ((fHS[i] + fBB[i])*r < -cutoffValue )? -cutoffValue : -(fHS[i] + fBB[i]);
            } else if ( type == 1) {
                printPotential = ( uBs1[i] > cutoffValue )? cutoffValue : uBs1[i];
                printForce = (fBs1[i]*r < -cutoffValue )? -cutoffValue : -fBs1[i];
            } else if ( type == 2) {
                printPotential = ( uBs2[i] > cutoffValue )? cutoffValue : uBs2[i];
                printForce = (fBs2[i]*r < -cutoffValue )? -cutoffValue : -fBs2[i];
            } else if ( type == 3) {
                printPotential = ( us1s2[i] > cutoffValue )? cutoffValue : us1s2[i];
                printForce = (fs1s2[i]*r < -cutoffValue )? -cutoffValue : -fs1s2[i];
            } else if ( type == 4) {
                printPotential = ( us1s1[i] > cutoffValue )? cutoffValue : us1s1[i];
                printForce = (fs1s1[i]*r < -cutoffValue )? -cutoffValue : -fs1s1[i];
            } else if ( type == 5) {
                printPotential = ( us2s2[i] > cutoffValue )? cutoffValue : us2s2[i];
                printForce = (fs2s2[i]*r < -cutoffValue )? -cutoffValue : -fs2s2[i];
            }
            // finally, you can print
            potentialOutputFile << printPotential << "\t" << printForce << "\n";
        }
        potentialOutputFile.close();
    }

}



int main( int argc, char *argv[] ) {
    // parse arguments
    if (argc != 3)
        error("Wrong number of arguments.\n");
    const std::string inputFileName(argv[1]);
    const std::string outputDirName(argv[2]);

    PotentialForLammps potential(inputFileName);
    potential.computePotential();
    potential.printPotentialsToFile(outputDirName);
}
