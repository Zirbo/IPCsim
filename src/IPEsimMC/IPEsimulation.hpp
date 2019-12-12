#include <fstream>
#include <vector>
#include "../IPC.hpp"
#include "../helpers/cell_lists.hpp"
#include "../helpers/randomNumberGenerator.hpp"
#include "../helpers/isotropic_pair_correlation_function.hpp"




struct SimulationStage {
    double inputTemperature;
    unsigned long inputStageTotalDuration;
    bool inputRestoringPreviousSimulation;
    bool inputPrintTrajectoryAndCorrelations;

    double deltaTrans;
    double deltaRot;

    SimulationStage() : inputTemperature{0.}, inputStageTotalDuration{0L},
                        inputRestoringPreviousSimulation{false}, inputPrintTrajectoryAndCorrelations{false},
                        deltaTrans{0.0}, deltaRot{0.0} {}
 };

class IPEsimulation {
public:
    IPEsimulation(SimulationStage const& stage);
    void run();


private:
    IPEsimulation();

    unsigned long simulationTime;
    std::ofstream outputFile;
    std::ofstream trajectoryFile;
    std::ofstream energyTrajectoryFile;

    bool printTrajectoryAndCorrelations;

    // state point
    int nIPEs;
    double density, temperature;
    // simulation duration
    double simulationTotalDuration;
    unsigned long printingInterval;
    // potential
    double e_BB, e_Bs, e_ss, e_min;
    double deltaOverSigma;
    double patchEccentricity;
    // work parameters
    double deltaPotential;
    double ipcRadius, ipcDiameter, patchRadius;
    double inverseTemperature;
    double ipcDiameterSquared;
    double potentialEnergy, simulationBoxSide;
    double BBinteractionRange, BBsquaredInteractionRange;
    double BsinteractionRange, BsSquaredInteractionRange;
    double ssInteractionRange, ssSquaredInteractionRange;
    double deltaTrans, deltaRot;
    double coeff_BB, coeff_Bs, coeff_ss;
    //double patchDistance, squaredPatchDistance, inversePatchDistance;
    // accounting
    double minimumSquaredDistance;

    std::vector<IPE> particles;
    cell_lists cells;

    IsotropicPairCorrelationFunction pairCorrelation;

    // selfexplanatory
    void initializeSystem(SimulationStage const& stage);
    void readInputFile();
    void printInputFileToOutputFile();
    void restorePreviousConfiguration();
    void initializeNewConfiguration();
    void computeSimulationStep();

    void makeRotationOrTranslationMove(IPE & ipe, RandomNumberGenerator &ranGen);
    double computeTotalPotential();
    void evaluateError(IPE const& firstIPE, IPE const& secndIPE);
    // the next functions return true if an overlap was detected, in which case dU is not to be used!!!
    const std::list<int> findAllTheIPEsInRange(IPE const& ipe);
    bool computePotentialOfAnIPEmove(IPE const& ipe, double& dU);
    bool computeInteractionsWithIPEsInList(const IPE &ipe, std::list<int> const& listOfIPEs, double& dU);
    bool computeInteractionsBetweenTwoIPEs(IPE const& firstIPE, IPE const& secndIPE, double& U);
    bool detectOverlap(const IPE &firstIPE, const IPE &secndIPE, const double rSquared);
    double computePotentialBetweenTwoIPEsInsideRange(const IPE &firstIPE, const IPE &secndIPE, const double r);
    const double computeOmega(const double Ra, const double Rb, const double rab);

    void outputSystemTrajectory(std::ofstream & outputTrajectoryFile);
    void outputSystemEnergies(std::ofstream &energyTrajectoryFile);

    // boundary conditions
    inline void absolutePBC(double & x) {  x -= std::floor(x);  }
    inline void relativePBC(double & x) {  x -= std::round(x);  }
    // Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
    void generateRandomOrientation(double (&a)[3], RandomNumberGenerator & r);
};
