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


#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/SpinorBasis/Spin.hpp"
#include "Basis/SpinorBasis/USpinorBasis.hpp"
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
     *  Create a spin-resolved ONV from textual/string representations.
     * 
     *  @param string_representation_alpha              the textual representation of the alpha-part of the spin-resolved ONV, for example "0011", indicating that the first two alpha-spin-orbitals should be occupied
     *  @param string_representation_beta               the textual representation of the beta-part of the spin-resolved ONV, for example "0011", indicating that the first two beta-spin-orbitals should be occupied
     * 
     *  @return a spin-unresolved ONV from textual/string representations.
     */
    static SpinResolvedONV FromString(const std::string& string_representation_alpha, const std::string& string_representation_beta);

    /**
     *  Create a spin-resolved ONV that represents the RHF single Slater determinant, occupying the N_P lowest alpha- and beta-spin-orbitals.
     * 
     *  @param K            the number of spatial orbitals
     *  @param N_P          the number of electron pairs
     * 
     *  @return a spin-resolved ONV that represents the RHF single Slater determinant
     * 
     * @note The ordering of the spin-orbitals is implicit: this method assumes that the spin-orbitals in the corresponding RSpinorBasis are sorted with increasing one-particle energy.
     */
    static SpinResolvedONV RHF(const size_t K, const size_t N_P);

    /**
     *  Create a spin-resolved ONV that represents the UHF single Slater determinant, occupying the N_alpha lowest alpha-spin-orbitals, and the N_beta lowest beta-spin-orbitals.
     * 
     *  @param K                the number of spatial orbitals
     *  @param N_alpha          the number of alpha-electrons
     *  @param N_beta           the number of beta-electrons
     * 
     *  @return a spin-resolved ONV that represents the UHF single Slater determinant
     * 
     * @note The ordering of the spin-orbitals is implicit: this method assumes that the spin-orbitals in the corresponding USpinorBasis are sorted with increasing one-particle energy.
     */
    static SpinResolvedONV UHF(const size_t K, const size_t N_alpha, const size_t N_beta);


    // PUBLIC METHODS

    /**
     *  Return a textual representation of this spin-resolved ONV.
     */
    std::string asString() const;

    /**
     *  Calculate the overlap <on|of>: the projection of between this spin-resolved ONV ('of') and another spin-resolved ONV ('on'), expressed in different R/U-spinor bases. The 'on'-ONV is supposed to be expressed in restricted spin-orbitals, and the 'of'-ONV is supposed to be expressed in unrestricted spin-orbitals.
     * 
     *  @param onv_on                       the spin-resolved ONV that should be projected on
     *  @param C_unrestricted               the coefficient matrix that describes the expansion of the alpha- and beta-spin-orbitals in terms of the underlying AOs
     *  @param C_restricted                 the coefficient matrix that describes the expansion of the restricted alpha/beta-spin-orbitals in terms of the underlying AOs
     *  @param S                            the overlap matrix of the underlying AOs
     * 
     *  @return the overlap element <on|of>
     * 
     *  @example This method can be used to project UHF-ONVs onto RHF-ONVs, by calling
     *          uhf_onv.calculateProjection(rhf_onv, C_unrestricted, C_restricted, S)
     */
    double calculateProjection(const SpinResolvedONV& onv_on, const SpinResolvedTransformationMatrix<double>& C_unrestricted, const TransformationMatrix<double>& C_restricted, const QCMatrix<double>& S) const;

    /**
     *  @param sigma                alpha or beta
     * 
     *  @return the number of sigma-electrons this spin-resolved ONV describes
     */
    size_t numberOfElectrons(const Spin sigma) const;

    /**
     *  @param sigma                alpha or beta
     * 
     *  @return the number of sigma-spatial orbitals/spin-orbitals that this ONV is expressed with
     */
    size_t numberOfSpatialOrbitals(const Spin sigma) const;

    /**
     *  @param sigma                alpha or beta
     * 
     *  @return the ONV that describes the occupations of the sigma-spin orbitals.
     */
    const SpinUnresolvedONV& onv(const Spin sigma) const;
};


}  // namespace GQCP
