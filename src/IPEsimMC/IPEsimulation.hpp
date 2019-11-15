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

    SimulationStage() : inputTemperature{0.}, inputStageTotalDuration{0L},
                        inputRestoringPreviousSimulation{false}, inputPrintTrajectoryAndCorrelations{false} {}
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
    // force and potential tables computation
    std::vector<double> uBB, uBs1, uBs2, us1s2, us1s1, us2s2;
    //void compilePotentialTables();
    double computeOmega(double Ra, double Rb, double rab);

    bool printTrajectoryAndCorrelations;

    // state point
    int nIPEs;
    double density, temperature;
    // simulation duration
    double simulationTotalDuration;
    unsigned long printingInterval;
    // potential
    double e_BB, e_Bs1, e_Bs2, e_s1s2, e_s1s1, e_s2s2, e_min;
    double ipcRadius, firstPatchRadius, firstPatchEccentricity, secndPatchRadius, secndPatchEccentricity;
    double forceAndEnergySamplingStep;
    // work parameters
    double potentialEnergy, simulationBoxSide;
    double interactionRange, squaredInteractionRange;
    //double patchDistance, squaredPatchDistance, inversePatchDistance;

    std::vector<IPE> particles;
    cell_lists cells;

    IsotropicPairCorrelationFunction pairCorrelation;

    // selfexplanatory
    void initializeSystem(SimulationStage const& stage);
    void readInputFile();
    void restorePreviousConfiguration();
    void initializeNewConfiguration();
    void computeSimulationStep();

    void makeRotationOrTranslationMove(IPE & ipe, RandomNumberGenerator &ranGen);
    double computePotentialDifference(IPE const& ipe);
    double computeInteractionsWithIPEsInTheSameCell(const IPE &ipe, std::list<int> const& ipesInCurrentCell);
    double computeInteractionsWithIPEsInNeighbouringCells(IPE const& ipe, std::list<int> const& ipesInNeighbouringCells);
    double computeInteractionsBetweenTwoIPEs(const int firstIPE, const int secndIPE);

    void computePotential();
    void outputSystemTrajectory(std::ofstream & outputTrajectoryFile);
    void outputSystemEnergies(std::ofstream &energyTrajectoryFile);

    // boundary conditions
    inline void absolutePBC(double & x) {  x -= std::floor(x);  }
    inline void relativePBC(double & x) {  x -= std::round(x);  }
    // Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
    void generateRandomOrientation(double (&a)[3], RandomNumberGenerator & r);
};


void test();
