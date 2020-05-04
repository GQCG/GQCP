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


#include "Basis/SpinorBasis/Spin.hpp"
#include "ONVBasis/SpinUnresolvedONV.hpp"


namespace GQCP {


/**
 *  An occupation number vector that is spin-resolved into alpha- and beta-constituents.
 */
class SpinResolvedONV {
private:
    SpinUnresolvedONV onv_alpha;  // the ONV that describes the occupations of the alpha spin-orbitals
    SpinUnresolvedONV onv_beta;   // the ONV that describes the occupations of the beta spin-orbitals

public:
    // CONSTRUCTORS

    /**
     *  @param onv_alpha                the ONV that describes the occupations of the alpha spin-orbitals
     *  @param onv_beta                 the ONV that describes the occupations of the beta spin-orbitals
     */
    SpinResolvedONV(const SpinUnresolvedONV& onv_alpha, const SpinUnresolvedONV& onv_beta);


    // OPERATORS

    /**
     *  @param other    the other spin-resolved ONV
     *
     *  @return if this spin-resolved ONV is the same as the other spin-resolved ONV
     */
    bool operator==(const SpinResolvedONV& other) const;

    /**
     *  @param other    the other spin-resolved ONV
     *
     *  @return if this spin-resolved ONV is not the same as the other spin-resolved ONV
     */
    bool operator!=(const SpinResolvedONV& other) const;


    // NAMED CONSTRUCTORS

    /**
     *  Create a spin-resolved ONV that represents the RHF single Slater determinant.
     * 
     *  @param K            the number of spatial orbitals
     *  @param N_P          the number of electron pairs
     * 
     *  @param a spin-resolved ONV that represents the RHF single Slater determinant
     */
    static SpinResolvedONV RHF(const size_t K, const size_t N_P);

    /**
     *  Create a spin-resolved ONV that represents the UHF single Slater determinant.
     * 
     *  @param K            the number of spatial orbitals
     *  @param N_P          the number of electron pairs
     * 
     *  @param a spin-resolved ONV that represents the UHF single Slater determinant
     */
    static SpinResolvedONV UHF(const size_t K, const size_t N_alpha, const size_t N_beta);


    // PUBLIC METHODS

    /**
     *  @param sigma                alpha or beta
     * 
     *  @return the ONV that describes the occupations of the sigma-spin orbitals.
     */
    const SpinUnresolvedONV& onv(Spin sigma) const;
};


}  // namespace GQCP
