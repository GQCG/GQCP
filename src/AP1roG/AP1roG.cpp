#include "AP1roG.hpp"


namespace GQCG {

/*
 *  CONSTRUCTORS
 */
/**
 *  Default constructor setting everything to zero
 */
AP1roG::AP1roG() :
    geminal_coefficients (AP1roGGeminalCoefficients()),
    electronic_energy (0.0)
{}


/**
 *  Constructor based on given @param geminal_coefficients and @param electronic_energy
 */
AP1roG::AP1roG(const GQCG::AP1roGGeminalCoefficients& geminal_coefficients, double electronic_energy) :
    geminal_coefficients (geminal_coefficients),
    electronic_energy (electronic_energy)
{}


}  // namespace GQCG
