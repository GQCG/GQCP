// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#include "Basis/ScalarBasis/GTOShell.hpp"

#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Utilities/miscellaneous.hpp"

#include <algorithm>


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  Create a `GTOShell` from all its specifications.
 * 
 *  @param l                                                        The angular momentum of the shell, i.e. the sum of the Cartesian exponents of this shell.
 *  @param nucleus                                                  The nucleus on which the shell is centered.
 *  @param gaussian_exponents                                       The Gaussian exponents for this shell, i.e. the exponents for the exponential. These are shared for every contraction.
 *  @param contraction_coefficients                                 The contraction coefficients for this shell.
 *  @param pure                                                     If this shell is considered to be 'pure', i.e. a spherical shell is pure, a Cartesian shell is not pure.
 *  @param are_embedded_normalization_factors_of_primitives         If the normalization factors of the primitives are embedded in the contraction coefficients.
 *  @param is_normalized                                            If the total normalization factor of the contracted GTO is already embedded in the contraction coefficients.
 */
GTOShell::GTOShell(const size_t l, const Nucleus& nucleus, const std::vector<double>& gaussian_exponents, const std::vector<double>& contraction_coefficients, const bool pure, const bool are_embedded_normalization_factors_of_primitives, const bool is_normalized) :
    pure {pure},
    are_embedded_normalization_factors_of_primitives {are_embedded_normalization_factors_of_primitives},
    normalized {is_normalized},
    l {l},
    m_nucleus {nucleus},
    gaussian_exponents {gaussian_exponents},
    contraction_coefficients {contraction_coefficients} {

    if (gaussian_exponents.size() != contraction_coefficients.size()) {
        throw std::invalid_argument("GTOShell(size_t, Nucleus, std::vector<double>, std::vector<double>): the exponents and contraction coefficients must match in size.");
    }
}


/*
 *  Shell characteristics
 */

/**
 *  @param rhs      the right-hand side of the operator ==
 *
 *  @return if this shell is considered equal to the other
 */
bool GTOShell::operator==(const GTOShell& rhs) const {

    /**
     *  A functor to compare two doubles with respect to a tolerance.
     */
    struct approx {
    public:
        double tolerance;


    public:
        approx(const double tolerance = 1.0e-12) :
            tolerance {tolerance} {}

        bool operator()(double lhs, double rhs) const {
            return std::abs(lhs - rhs) < tolerance;
        }
    };


    // Return if all members are equal/close
    return (this->l == rhs.l) &&
           (Nucleus::equalityComparer()(this->m_nucleus, rhs.m_nucleus)) &&
           (this->pure == rhs.pure) &&
           (std::equal(this->gaussian_exponents.begin(), this->gaussian_exponents.end(), rhs.gaussian_exponents.begin(), approx())) &&
           (std::equal(this->contraction_coefficients.begin(), this->contraction_coefficients.end(), rhs.contraction_coefficients.begin(), approx()));
}


/*
 *  MARK: Normalization
 */

/**
 *  Embed the total normalization factor of the corresponding linear combination of spherical (or axis-aligned Cartesian) GTOs into the contraction coefficients.
 */
void GTOShell::embedNormalizationFactor() {

    if (!this->normalized) {

        // Calculate the total norm of the shell's basis functions (corresponding to axis-aligned Cartesian GTOs).
        double norm = 0.0;
        for (size_t i = 0; i < this->contractionSize(); i++) {
            double c_i = this->contraction_coefficients[i];
            double alpha_i = this->gaussian_exponents[i];

            for (size_t j = 0; j < this->contractionSize(); j++) {
                double c_j = this->contraction_coefficients[j];
                double alpha_j = this->gaussian_exponents[j];

                double pair_exponent = (alpha_i + alpha_j) / 2;
                norm += c_i * c_j * CartesianGTO::calculateNormalizationFactor(pair_exponent, CartesianExponents(this->l, 0, 0));
            }
        }

        // Multiply all the contraction coefficients with the total normalization factor.
        for (size_t i = 0; i < this->contractionSize(); i++) {
            this->contraction_coefficients[i] *= std::pow(norm, -1.0 / 2.0);
        }
        this->normalized = true;
    }
}


/**
 *  Embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing.
 *
 *  @note The normalization factor that is embedded, corresponds to the spherical (or axis-aligned Cartesian) GTO.
 */
void GTOShell::embedNormalizationFactorsOfPrimitives() {

    if (!this->are_embedded_normalization_factors_of_primitives) {
        for (size_t i = 0; i < this->contractionSize(); i++) {
            this->contraction_coefficients[i] *= CartesianGTO::calculateNormalizationFactor(this->gaussian_exponents[i], CartesianExponents(this->l, 0, 0));  // The normalization factor of an axis-aligned Cartesian GTO.
        }

        this->are_embedded_normalization_factors_of_primitives = true;
    }
}


/**
 *  Un-embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing.
 *
 *  @note The the normalization factor that is un-embedded corresponds to the spherical (or axis-aligned Cartesian) GTO.
 */
void GTOShell::unEmbedNormalizationFactorsOfPrimitives() {

    if (this->are_embedded_normalization_factors_of_primitives) {
        for (size_t i = 0; i < this->contractionSize(); i++) {
            this->contraction_coefficients[i] /= CartesianGTO::calculateNormalizationFactor(this->gaussian_exponents[i], CartesianExponents(this->l, 0, 0));  // The normalization factor of an axis-aligned Cartesian GTO.
        }

        this->are_embedded_normalization_factors_of_primitives = false;
    }
}


/*
 *  MARK: Basis functions
 */

/**
 *  @return The number of basis functions that this shell contains.
 */
size_t GTOShell::numberOfBasisFunctions() const {

    if (pure) {  // If the shell is spherical.
        return 2 * this->l + 1;
    } else {  // If the shell is Cartesian.
        return (this->l + 1) * (this->l + 2) / 2;
    }
}


/**
 *  @return A list of the Cartesian exponents that have this shell's angular momentum (in lexicographical ordering).
 */
std::vector<CartesianExponents> GTOShell::generateCartesianExponents() const {

    // The different Cartesian exponents are the ways the angular momentum can be divided in 3 (x,y,z) partitions.
    const auto partitions = generatePartitionsOf(this->angularMomentum(), 3);


    // Cast the triples into CartesianExponents.
    std::vector<CartesianExponents> all_cartesian_exponents;
    all_cartesian_exponents.reserve(partitions.size());

    std::transform(partitions.begin(), partitions.end(),
                   std::back_inserter(all_cartesian_exponents),
                   [](const std::vector<size_t>& partition) {
                       return CartesianExponents(partition);
                   });

    return all_cartesian_exponents;
}


/**
 *  Construct all basis functions contained in this shell.
 * 
 *  @return The basis functions that correspond to this shell.
 * 
 *  @note The basis functions are ordered lexicographically. This means x < y < z.
 */
std::vector<GTOShell::BasisFunction> GTOShell::basisFunctions() const {

    // Prepare some variables.
    const auto& contraction_coefficients = this->contractionCoefficients();
    const auto& gaussian_exponents = this->gaussianExponents();


    // Generate the Cartesian exponents in a lexicographical ordering.
    const auto all_cartesian_exponents = this->generateCartesianExponents();


    // Do the actual 'contraction' of the primitives and the contraction coefficients.
    std::vector<GTOShell::BasisFunction> basis_functions;
    basis_functions.reserve(this->numberOfBasisFunctions());
    for (const auto& cartesian_exponents : all_cartesian_exponents) {

        GTOShell::BasisFunction basis_function;
        for (size_t d = 0; d < this->contractionSize(); d++) {
            auto coefficient = contraction_coefficients[d];
            const CartesianGTO function {gaussian_exponents[d], cartesian_exponents, this->nucleus().position()};

            basis_function.append({coefficient}, {function});
        }

        basis_functions.push_back(basis_function);
    }


    return basis_functions;
}


}  // namespace GQCP
