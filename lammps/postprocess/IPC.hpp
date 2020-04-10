#ifndef __IPC_HEADER_INCLUDED__
#define __IPC_HEADER_INCLUDED__

#include <vector>
#include <array>

struct Particle {
    double x[3], v[3], F[3];
};

struct IPCbase {
  Particle ipcCenter;
  char type;
  int number;
};

struct IPC : public IPCbase {
    Particle firstPatch, secndPatch;
    double eFp1[3], eFp2[3];
};

typedef std::vector<IPC> Ensemble;
typedef std::vector<std::array<double, 3>> SpaceVector;

#endif // __IPC_HEADER_INCLUDED__
