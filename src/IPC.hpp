#ifndef __IPC_HEADER_INCLUDED__
#define __IPC_HEADER_INCLUDED__

struct Particle {
    double x[3], v[3], F[3];
    // TODO: when all the sim is working, remove the forces from here and only keep the IPC effective forces
};
struct IPC {
    Particle ipcCenter, firstPatch, secndPatch;
    double eFp1[3], eFp2[3];
    char type;
    // unfortunately I need the number for managing the lists...
    // there's probably another way but I don't see it now, and anyway this is cheap
    int number;
};

#endif // __IPC_HEADER_INCLUDED__
