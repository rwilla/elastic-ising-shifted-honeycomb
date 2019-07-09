#include "lattice.hpp"

namespace lattice {

NearestNeighbors::NearestNeighbors(int N){
    this->AB1 = std::vector<int>(N);
    this->AB2 = std::vector<int>(N);
    this->AB3 = std::vector<int>(N);
    this->AB4 = std::vector<int>(N);
    this->AB5 = std::vector<int>(N);
    this->AB6 = std::vector<int>(N);
    this->AA1 = std::vector<int>(N);
    this->AA2 = std::vector<int>(N);
    this->BA1 = std::vector<int>(N);
    this->BA2 = std::vector<int>(N);
    this->BA3 = std::vector<int>(N);
    this->BA4 = std::vector<int>(N);
    this->BA5 = std::vector<int>(N);
    this->BA6 = std::vector<int>(N);
    this->BB1 = std::vector<int>(N);
    this->BB2 = std::vector<int>(N);
    this->sweeporder = std::vector<int>(N);
}
    
}
