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


#include "Basis/TransformationMatrix.hpp"
#include "ONVBasis/BaseONVBasis.hpp"
#include "Mathematical/Representation/Matrix.hpp"


namespace GQCP {


/**
 *  A class that represents a wave function: expansion coefficients in a ONV basis
 */
class LinearExpansion {
private:
    std::shared_ptr<BaseONVBasis> fock_space;
    VectorX<double> coefficients;  // the expansion coefficients in the ONV basis

public:
    // CONSTRUCTORS
    LinearExpansion() = default;

    /**
     *  Construct a normalized wave function from possibly non-normalized coefficients
     *
     *  @param base_fock_space      the ONV basis in which the wave function 'lives'
     *  @param coefficients         the expansion coefficients
     */
    LinearExpansion(const BaseONVBasis& base_fock_space, const VectorX<double>& coefficients);


    // GETTERS
    const VectorX<double>& get_coefficients() const { return coefficients; }
    const BaseONVBasis& get_fock_space() const { return *fock_space; }


    // PUBLIC METHODS
    /**
     *  @return the Shannon entropy (or information content) of the wave function
     */
    double calculateShannonEntropy() const;

    /**
     *  Transform the underlying ONV basis of the wave function (only for FCI [ProductONVBasis]) and recalculate the ONV expansion coefficients
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
     void basisTransform(const TransformationMatrix<double>& T);

    /** 
     *  @param other            wave function for the comparison
     *  @param tolerance        tolerance for the comparison of coefficients
     * 
     *  @return if two wave functions are equal within a given tolerance
     */
     bool isApprox(const LinearExpansion& other, double tolerance = 1e-10) const;
};


}  // namespace GQCP
