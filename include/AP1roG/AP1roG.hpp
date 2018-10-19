// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#ifndef AP1roG_hpp
#define AP1roG_hpp


#include "AP1roGGeminalCoefficients.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"


namespace GQCP {


/**
 *  A class representing an AP1roG wave function: it holds the geminal coefficients that are a solution to the AP1roG projected Schrödinger equations
 */
class AP1roG {
private:
    GQCP::AP1roGGeminalCoefficients geminal_coefficients;

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
    AP1roG(const GQCP::AP1roGGeminalCoefficients& geminal_coefficients, double electronic_energy);


    // GETTERS
    GQCP::AP1roGGeminalCoefficients get_geminal_coefficients() const { return this->geminal_coefficients; }
    double get_electronic_energy() const { return this->electronic_energy; }
};


/*
 *  HELPER FUNCTIONS
 */
/**
 *  Calculate the AP1roG energy given AP1roG geminal coefficients @param G and Hamiltonian parameters @param ham_par
 */
double calculateAP1roGEnergy(const GQCP::AP1roGGeminalCoefficients& G, const GQCP::HamiltonianParameters& ham_par);


}  // namespace GQCP


#endif /* AP1roG_hpp */
