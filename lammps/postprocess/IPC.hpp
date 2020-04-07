#ifndef __IPC_HEADER_INCLUDED__
#define __IPC_HEADER_INCLUDED__

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


#endif // __IPC_HEADER_INCLUDED__
