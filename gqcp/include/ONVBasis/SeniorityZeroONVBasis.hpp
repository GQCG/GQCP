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


#include "Mathematical/Representation/Matrix.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"


namespace GQCP {


/**
 *  A spin-resolved ONV basis that contains all seniority-zero (i.e. doubly-occupied) (spin-resolved) ONVs.
 */
class SeniorityZeroONVBasis {
private:
    size_t K;  // the number of spatial orbitals
    size_t N_P;  // the number of electron pairs

    size_t dim; // the dimension of this ONV basis


public:

    // CONSTRUCTORS

    /**
     *  @param K            the number of spatial orbitals
     *  @param N_P          the number of electron pairs
     */
    SeniorityZeroONVBasis(const size_t K, const size_t N_P);


    // STATIC PUBLIC METHODS

    /**
     *  @param K            the number of spatial orbitals
     *  @param N_P          the number of electron pairs
     * 
     *  @return the dimension of a seniority-zero ONV basis with the given number of spatial orbitals and electron pairs
     */
    static size_t calculateDimension(const size_t K, const size_t N_P);


    // PUBLIC METHODS

    /**
     *  @return the dimension of this ONV basis
     */
    size_t dimension() const { return this->dim; }

    /**
     *  @return a coefficient vector that describes the expansion coefficients of the Hartree-Fock wave function (i.e. the single Slater determinant with the lowest energy)
     */
    VectorX<double> hartreeFockExpansion() const;

    /**
     *  @return the number of electron pairs that this ONV basis is related to
     */
    size_t numberOfElectronPairs() const { return this->N_P; }

    /**
     *  @return the number of spatial orbitals that this ONV basis is related to
     */
    size_t numberOfSpatialOrbitals() const { return this->K; }

    /**
     *  @return a spin-unresolved ONV basis that behaves analogously (with respect to a doubly-occupied situation) as this seniority-zero ONV basis
     */
    SpinUnresolvedONVBasis proxy() const { return SpinUnresolvedONVBasis(this->K, this->N_P); }
};


}  // namespace GQCP
