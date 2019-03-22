// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
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
#include "Basis/BasisFunction.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param lc       a linear combination of CartesianGTOs
 */
BasisFunction::BasisFunction(const Base& lc) :
    Base(lc)
{

    // Check if all the CartesianGTOs have the same Cartesian exponents
    auto ref_cartesian_exponents = lc.get_functions()[0].get_cartesian_exponents();
    for (const auto& gto : lc.get_functions()) {
        if (gto.get_cartesian_exponents() != ref_cartesian_exponents) {
            throw std::invalid_argument("BasisFunction::BasisFunction(Base&): all the CartesianGTOs should have the same Cartesian exponents");
        }
    }

    // Modify the contraction coefficients (i.e. the coefficients in the linear combination) by the total normalization factor in order to obtain a normalized contracted GTO
    this->N_total = this->calculateNormalizationFactor();
    for (auto& contraction_coefficient : this->coefficients) {
        contraction_coefficient *= this->N_total;
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the total normalization factor of this basis function
 */
double BasisFunction::calculateNormalizationFactor() const {

    // KISS-implementation of the total normalization factor

    double sum = 0.0;
    for (size_t a = 0; a < this->length(); a++) {
        double value = 1.0;

        double d_a = this->coefficients[a];
        double N_a = this->functions[a].get_N();
        double alpha_a = this->functions[a].get_gaussian_exponent();

        for (size_t b = 0; b < this->length(); b++) {
            double d_b = this->coefficients[b];
            double N_b = this->functions[b].get_N();
            double alpha_b = this->functions[b].get_gaussian_exponent();

            value *= N_a * d_a * N_b * d_b;

            for (const auto& direction : {CartesianDirection::x, CartesianDirection::y, CartesianDirection::z}) {

                size_t cartesian_exponent = this->functions[a].get_cartesian_exponents().value(direction);  // same for every primitive

                value *= CartesianGTO::calculateNormalizationFactorComponent((alpha_a + alpha_b)/2, cartesian_exponent);
            }

            sum += value;
        }
    }


    return std::pow(sum, -0.5);
}


}  // namespace GQCP
