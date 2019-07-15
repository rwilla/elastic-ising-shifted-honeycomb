#include "hamiltonian.hpp"
using namespace std;
namespace hamiltonian {

//EISH = elastic ising shifted honeycomb
namespace EISH {
// onsite
double onsite(const iparameters::InitialParameters& ip, double S, double D, int spin_index)
{
    double e = 0.0;
    
    //displacement energy
    e += ip.eps / 2 * D * D;
    //field-dependence
    e -= ip.Hz * S;
    
    return e;
}
//Full hamiltonian for a spin on sublattice A
double siteA(const iparameters::InitialParameters& ip, double S, double D, int spin_index, const lattice::EISH::NearestNeighbors& nn, observables::EISH::SystemState& s)
{
    double e = 0.0;
    e += hamiltonian::EISH::onsite(ip, S, D, spin_index);
    e += hamiltonian::EISH::nnA(ip, S, D, spin_index, nn, s);
    
    return e;
}
//Full hamiltonian for a spin on sublattice B
double siteB(const iparameters::InitialParameters& ip, double S, double D, int spin_index, const lattice::EISH::NearestNeighbors& nn, observables::EISH::SystemState& s)
{
    double e = 0.0;
    e += hamiltonian::EISH::onsite(ip, S, D, spin_index);
    e += hamiltonian::EISH::nnB(ip, S, D, spin_index, nn, s);
    
    return e;
}
//Nearest neighbor interaction for a spin on sublattice A
double nnA(const iparameters::InitialParameters& ip, double S, double D, int spin_index, const lattice::EISH::NearestNeighbors& nn, observables::EISH::SystemState& s)
{
    double e = 0.0;
        
    double nSB;
    double nDB;
            
    nSB = s.SB[nn.AB1[spin_index]];
    nDB = s.DB[nn.AB1[spin_index]];
    
    e += - ip.Jab * (1.0 + D - nDB) * (S * nSB);
    
    nSB = s.SB[nn.AB2[spin_index]];
    nDB = s.DB[nn.AB2[spin_index]];
    
    e += - ip.Jab * (1.0 + D - nDB) * (S * nSB);
    
    nSB = s.SB[nn.AB3[spin_index]];
    nDB = s.DB[nn.AB3[spin_index]];
    
    e += - ip.Jab * (1.0 + D - nDB) * (S * nSB);
    
    nSB = s.SB[nn.AB4[spin_index]];
    nDB = s.DB[nn.AB4[spin_index]];
    
    e += - ip.Jab * (1.0 - D + nDB) * (S * nSB);
    
    nSB = s.SB[nn.AB5[spin_index]];
    nDB = s.DB[nn.AB5[spin_index]];
    
    e += - ip.Jab * (1.0 - D + nDB) * (S * nSB);
    
    nSB = s.SB[nn.AB6[spin_index]];
    nDB = s.DB[nn.AB6[spin_index]];
    
    e += - ip.Jab * (1.0 - D + nDB) * (S * nSB);
    
    double nSA;
    double nDA;
    nSA = s.SA[nn.AA1[spin_index]];
    nDA = s.DA[nn.AA1[spin_index]];
    
    e += - ip.Jcc * (1.0 + D - nDA) * (S * nSA);
    
    nSA = s.SA[nn.AA2[spin_index]];
    nDA = s.DA[nn.AA2[spin_index]];
    
    e += - ip.Jcc * (1.0 - D + nDA) * (S * nSA);
    
    return e;
}

//Nearest neighbor interaction for a spin on sublattice B
double nnB(const iparameters::InitialParameters& ip, double S, double D, int spin_index, const lattice::EISH::NearestNeighbors& nn, observables::EISH::SystemState& s)
{
    double e = 0.0;
        
    double nSA;
    double nDA;
            
    nSA = s.SA[nn.BA1[spin_index]];
    nDA = s.DA[nn.BA1[spin_index]];
    
    e += - ip.Jab * (1.0 + D - nDA) * (S * nSA);
    
    nSA = s.SA[nn.BA2[spin_index]];
    nDA = s.DA[nn.BA2[spin_index]];
    
    e += - ip.Jab * (1.0 + D - nDA) * (S * nSA);
    
    nSA = s.SA[nn.BA3[spin_index]];
    nDA = s.DA[nn.BA3[spin_index]];
    
    e += - ip.Jab * (1.0 + D - nDA) * (S * nSA);
    
    nSA = s.SA[nn.BA4[spin_index]];
    nDA = s.DA[nn.BA4[spin_index]];
    
    e += - ip.Jab * (1.0 - D + nDA) * (S * nSA);
    
    nSA = s.SA[nn.BA5[spin_index]];
    nDA = s.DA[nn.BA5[spin_index]];
    
    e += - ip.Jab * (1.0 - D + nDA) * (S * nSA);
    
    nSA = s.SA[nn.BA6[spin_index]];
    nDA = s.DA[nn.BA6[spin_index]];
    
    e += - ip.Jab * (1.0 - D + nDA) * (S * nSA);
    
    double nSB;
    double nDB;
    nSB = s.SB[nn.BB1[spin_index]];
    nDB = s.DB[nn.BB1[spin_index]];
    
    e += - ip.Jcc * (1.0 + D - nDB) * (S * nSB);
    
    nSB = s.SB[nn.BB2[spin_index]];
    nDB = s.DB[nn.BB2[spin_index]];
    
    e += - ip.Jcc * (1.0 - D + nDB) * (S * nSB);
    
    return e;
}
}
}
