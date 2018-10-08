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


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on given @param geminal_coefficients
     */
    AP1roG(const GQCG::AP1roGGeminalCoefficients& geminal_coefficients);


    // GETTERS
    GQCG::AP1roGGeminalCoefficients get_geminal_coefficients() const { return this->geminal_coefficients; }
};



}  // namespace GQCG


#endif /* AP1roG_hpp */
