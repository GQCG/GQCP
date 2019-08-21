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
#pragma once


#include "Basis/CartesianGTO.hpp"
#include "Molecule/Nucleus.hpp"


namespace GQCP {


/**
 *  A class that represents a shell of GTOs: it specifies in a condensed way which GTOs are on an nucleus
 */
class GTOShell {
private:
    bool pure;  // true if spherical, false if Cartesian
    bool embedded_normalization_factors_of_primitives;  // if the normalization factors of the primitives are embedded in the contraction coefficients
    bool normalized;  // if the total normalization factor is already embedded in the contraction coefficients
    size_t l;  // the angular momentum of the shell
    Nucleus nucleus;  // nucleus on which the shell is centered
    std::vector<double> gaussian_exponents;  // Gaussian exponents (i.e. for the exponential), shared for every contraction
    std::vector<double> contraction_coefficients;


public:
    using BasisFunction = CartesianGTO;


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
    GTOShell(size_t l, const Nucleus& nucleus, const std::vector<double>& gaussian_exponents, const std::vector<double>& contraction_coefficients, bool pure=true, bool are_embedded_normalization_factors_of_primitives=false, bool is_normalized=false);


    // GETTERS
    bool is_pure() const { return this->pure; }
    bool are_embedded_normalization_factors_of_primitives() const { return this->embedded_normalization_factors_of_primitives; }
    bool is_normalized() const { return this->normalized; }
    size_t get_l() const { return this->l; }
    const Nucleus& get_nucleus() const { return this->nucleus; }
    const std::vector<double>& get_gaussian_exponents() const { return this->gaussian_exponents; }
    const std::vector<double>& get_contraction_coefficients() const { return this->contraction_coefficients; }


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
     *  @return the number of basis functions that are in this shell
     */
    size_t numberOfBasisFunctions() const;

    /**
     *  @return the size of the contraction in the shell, i.e. the number of primitives contracted in this shell
     */
    size_t contractionSize() const;

    /**
     *  Embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing
     *
     *  Note that the normalization factor that is embedded corresponds to the spherical (or axis-aligned Cartesian) GTO
     */
    void embedNormalizationFactorsOfPrimitives();

    /**
     *  Embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing
     *
     *  Note that the normalization factor that is embedded corresponds to the spherical (or axis-aligned Cartesian) GTO
     */
    void unEmbedNormalizationFactorsOfPrimitives();

    /**
     *  Embed the total normalization factor of the corresponding linear combination of spherical (or axis-aligned Cartesian) GTOs into the contraction coefficients
     */
    void embedNormalizationFactor();
};


}  // namespace GQCP
