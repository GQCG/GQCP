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

#pragma once


#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/LinearCombination.hpp"
#include "Molecule/Nucleus.hpp"


namespace GQCP {


/**
 *  A class that represents a shell of GTOs: it specifies in a condensed way which GTOs are on an nucleus.
 */
class GTOShell {
private:
    bool pure;                                          // true if spherical, false if Cartesian
    bool embedded_normalization_factors_of_primitives;  // if the normalization factors of the primitives are embedded in the contraction coefficients
    bool normalized;                                    // if the total normalization factor is already embedded in the contraction coefficients
    size_t l;                                           // the angular momentum of the shell
    Nucleus m_nucleus;                                  // nucleus on which the shell is centered
    std::vector<double> gaussian_exponents;             // Gaussian exponents (i.e. for the exponential), shared for every contraction
    std::vector<double> contraction_coefficients;


public:
    using Primitive = CartesianGTO;                              // the type of primitives that this shell is made up with
    using BasisFunction = LinearCombination<double, Primitive>;  // the type of basis functions that this shell can produce


public:
    // CONSTRUCTORS
    /**
     *  @param l                                                    the angular momentum of the shell
     *  @param nucleus                                              the nucleus on which the shell is centered
     *  @param gaussian_exponents                                   the Gaussian exponents, which are shared for every contraction
     *  @param contraction_coefficients                             the contraction coefficients
     *  @param pure                                                 whether the shell is considered to be spherical or not
     *  @param are_embedded_normalization_factors_of_primitives     if the normalization factors of the primitives are embedded in the contraction coefficients
     *  @param is_normalized                                        if the total normalization factor is already embedded in the contraction coefficients
     */
    GTOShell(const size_t l, const Nucleus& nucleus, const std::vector<double>& gaussian_exponents, const std::vector<double>& contraction_coefficients, const bool pure = true, const bool are_embedded_normalization_factors_of_primitives = false, const bool is_normalized = false);


    // OPERATORS
    /**
     *  @param rhs      the right-hand side of the operator ==
     *
     *  @return if this shell is considered equal to the other
     */
    bool operator==(const GTOShell& rhs) const;


    // PUBLIC METHODS

    /**
     *  @return the angular momentum of the shell
     */
    size_t angularMomentum() const { return this->l; }

    /**
     *  @return if the normalization factors of the primitives are embedded in the contraction coefficients
     */
    bool areEmbeddedNormalizationFactorsOfPrimitives() const { return this->embedded_normalization_factors_of_primitives; }

    /**
     *  @return the basis functions that correspond to this shell
     * 
     *  @note The basis functions are ordered lexicographically. This means x < y < z.
     */
    std::vector<BasisFunction> basisFunctions() const;

    /**
     *  @return the contraction coefficients for this shell
     */
    const std::vector<double>& contractionCoefficients() const { return this->contraction_coefficients; }

    /**
     *  @return the size of the contraction in the shell, i.e. the number of primitives contracted in this shell
     */
    size_t contractionSize() const { return this->contraction_coefficients.size(); }

    /**
     *  Embed the total normalization factor of the corresponding linear combination of spherical (or axis-aligned Cartesian) GTOs into the contraction coefficients.
     */
    void embedNormalizationFactor();

    /**
     *  Embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing.
     *
     *  @note The normalization factor that is embedded, corresponds to the spherical (or axis-aligned Cartesian) GTO
     */
    void embedNormalizationFactorsOfPrimitives();

    /**
     *  @return the Gaussian exponents (i.e. for the exponential) for this shell, shared for every contraction
     */
    const std::vector<double>& gaussianExponents() const { return this->gaussian_exponents; }

    /**
     *  @return a list of the Cartesian exponents that have this shell's angular momentum (in lexicographical ordering).
     */
    std::vector<CartesianExponents> generateCartesianExponents() const;

    /**
     *  @return if the total normalization factor is already embedded in the contraction coefficients
     */
    bool isNormalized() const { return this->normalized; }

    /**
     *  @return true if the GTO shell is spherical, false if it is Cartesian
     */
    bool isPure() const { return this->pure; }

    /**
     *  @return The nucleus on which the shell is centered.
     */
    const Nucleus& nucleus() const { return this->m_nucleus; }

    /**
     *  @return the number of basis functions that are in this shell
     */
    size_t numberOfBasisFunctions() const;

    /**
     *  Embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing.
     *
     *  @note The the normalization factor that is embedded corresponds to the spherical (or axis-aligned Cartesian) GTO
     */
    void unEmbedNormalizationFactorsOfPrimitives();
};


}  // namespace GQCP
