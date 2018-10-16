#ifndef AP1roG_hpp
#define AP1roG_hpp


#include "AP1roGGeminalCoefficients.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"


namespace GQCG {


/**
 *  A class representing an AP1roG wave function: it holds the geminal coefficients that are a solution to the AP1roG projected Schrödinger equations
 */
class AP1roG {
private:
    GQCG::AP1roGGeminalCoefficients geminal_coefficients;

    double electronic_energy;

public:
    // CONSTRUCTORS
    /**
     *  Default constructor setting everything to zero
     */
    AP1roG();

    /**
     *  Constructor based on given @param geminal_coefficients and @param electronic_energy
     */
    AP1roG(const GQCG::AP1roGGeminalCoefficients& geminal_coefficients, double electronic_energy);


    // GETTERS
    GQCG::AP1roGGeminalCoefficients get_geminal_coefficients() const { return this->geminal_coefficients; }
    double get_electronic_energy() const { return this->electronic_energy; }
};


/*
 *  HELPER FUNCTIONS
 */
/**
 *  Calculate the AP1roG energy given AP1roG geminal coefficients @param G and Hamiltonian parameters @param ham_par
 */
double calculateAP1roGEnergy(const GQCG::AP1roGGeminalCoefficients& G, const GQCG::HamiltonianParameters& ham_par);


}  // namespace GQCG


#endif /* AP1roG_hpp */
