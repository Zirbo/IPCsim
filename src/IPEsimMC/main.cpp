#include <iostream>
#include "IPEsimulation.hpp"

int main() {
    SimulationStage stage;
    stage.inputTemperature = 1.0;
    stage.inputStageTotalDuration = 1000;
    stage.inputPrintTrajectoryAndCorrelations = true;

    IPEsimulation simulation(stage);
    simulation.run();
}
