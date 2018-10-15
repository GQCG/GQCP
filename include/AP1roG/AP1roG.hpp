#ifndef AP1roG_hpp
#define AP1roG_hpp


#include "AP1roGGeminalCoefficients.hpp"


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

    /**
     *  Copy assignment
     */
    AP1roG& operator=(const AP1roG& ap1rog) = default;


    // GETTERS
    GQCG::AP1roGGeminalCoefficients get_geminal_coefficients() const { return this->geminal_coefficients; }
    double get_electronic_energy() const { return this->electronic_energy; }
};



}  // namespace GQCG


#endif /* AP1roG_hpp */
