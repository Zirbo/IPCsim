#ifndef __IPC_HEADER_INCLUDED__
#define __IPC_HEADER_INCLUDED__

struct Particle {
    double x[3], v[3], F[3];
};
struct IPC {
    Particle center, firstPatch, secndPatch;
    double eFp1[3], eFp2[3];
    char type;
};

#endif // __IPC_HEADER_INCLUDED__
