#ifndef __IPCPOSTPROCESS_POTENTIAL_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_POTENTIAL_HEADER_INCLUDED__

#include <string>
#include <vector>

class IPCpotential {
public:
    double spacing;
    std::vector<double> uBB, uBs, uss;
    void initialize(std::string const& potentialsDir);

private:
    void readFile(std::string const& type, double &spacing, std::vector<double> & container);
};

#endif //__IPCPOSTPROCESS_POTENTIAL_HEADER_INCLUDED__
