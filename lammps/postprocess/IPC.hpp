#ifndef __IPC_HEADER_INCLUDED__
#define __IPC_HEADER_INCLUDED__

#include <vector>
#include <array>

struct Particle {
    double x[3], v[3], F[3];
};

struct IPC {
  Particle ipcCenter;
  Particle firstPatch, secndPatch;
  double eFp1[3], eFp2[3];
  char type;
  int number;
};

typedef std::vector<IPC> Ensemble;
typedef std::array<double, 3> Triad;
typedef std::vector<std::array<double, 3>> VectorOfTriads;

const Triad DIMENSIONS = {0, 1, 2};

#endif // __IPC_HEADER_INCLUDED__
