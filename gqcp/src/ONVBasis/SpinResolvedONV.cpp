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

#include "ONVBasis/SpinResolvedONV.hpp"

#include "ONVBasis/SpinUnresolvedONVBasis.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param onv_alpha                the ONV that describes the occupations of the alpha spin-orbitals
 *  @param onv_beta                 the ONV that describes the occupations of the beta spin-orbitals
 */
SpinResolvedONV::SpinResolvedONV(const SpinUnresolvedONV& onv_alpha, const SpinUnresolvedONV& onv_beta) :
    onv_alpha {onv_alpha},
    onv_beta {onv_beta} {}


/*
 *  OPERATORS
 */

/**
 *  @param other    the other spin-resolved ONV
 *
 *  @return if this spin-resolved ONV is the same as the other spin-resolved ONV
 */
bool SpinResolvedONV::operator==(const SpinResolvedONV& other) const {

    return (this->onv(Spin::alpha) == other.onv(Spin::alpha)) && (this->onv(Spin::beta) == other.onv(Spin::beta));
}


/**
 *  @param other    the other spin-resolved ONV
 *
 *  @return if this spin-resolved ONV is not the same as the other spin-resolved ONV
 */
bool SpinResolvedONV::operator!=(const SpinResolvedONV& other) const {

    return !(*this == other);
}


/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  Create a spin-resolved ONV from textual/string representations.
 * 
 *  @param string_representation_alpha              the textual representation of the alpha-part of the spin-resolved ONV, for example "0011", indicating that the first two alpha-spin-orbitals should be occupied
 *  @param string_representation_beta               the textual representation of the beta-part of the spin-resolved ONV, for example "0011", indicating that the first two beta-spin-orbitals should be occupied
 * 
 *  @return a spin-unresolved ONV from textual/string representations.
 */
SpinResolvedONV SpinResolvedONV::FromString(const std::string& string_representation_alpha, const std::string& string_representation_beta) {

    const auto onv_alpha = SpinUnresolvedONV::FromString(string_representation_alpha);
    const auto onv_beta = SpinUnresolvedONV::FromString(string_representation_beta);

    return SpinResolvedONV(onv_alpha, onv_beta);
}


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
SpinResolvedONV SpinResolvedONV::RHF(const size_t K, const size_t N_P) {

    const auto sigma_onv = SpinUnresolvedONV::GHF(K, N_P);  // for sigma = alpha and sigma = beta

    return SpinResolvedONV(sigma_onv, sigma_onv);
}


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
SpinResolvedONV SpinResolvedONV::UHF(const size_t K, const size_t N_alpha, const size_t N_beta) {

    const auto alpha_onv = SpinUnresolvedONV::GHF(K, N_alpha);
    const auto beta_onv = SpinUnresolvedONV::GHF(K, N_beta);

    return SpinResolvedONV(alpha_onv, beta_onv);
}


/*
 *  PUBLIC METHODS
 */


/**
 *  @return a textual representation of this spin-resolved ONV.
 */
std::string SpinResolvedONV::asString() const {

    return this->onv(Spin::alpha).asString() + "|" + this->onv(Spin::beta).asString();
}


/**
 *  Calculate the overlap <on|of>: the projection of between this spin-resolved ONV ('of') and another spin-resolved ONV ('on'), expressed in different R/U-spinor bases. The 'on'-ONV is supposed to be expressed in restricted spin-orbitals, and the 'of'-ONV is supposed to be expressed in unrestricted spin-orbitals.
 * 
 *  @param onv_on                       the spin-resolved ONV that should be projected on
 *  @param C_alpha                      the coefficient matrix that describes the expansion of the alpha-spin-orbitals in terms of the underlying AOs
 *  @param C_beta                       the coefficient matrix that describes the expansion of the beta-spin-orbitals in terms of the underlying AOs
 *  @param C                            the coefficient matrix that describes the expansion of the restricted alpha/beta-spin-orbitals in terms of the underlying AOs
 *  @param S                            the overlap matrix of the underlying AOs
 * 
 *  @return the overlap element <on|of>
 * 
 *  @example This method can be used to project UHF-ONVs onto RHF-ONVs, by calling
 *          uhf_onv.calculateProjection(rhf_onv, C_alpha, C_beta, C, S)
 */
double SpinResolvedONV::calculateProjection(const SpinResolvedONV& onv_on, const TransformationMatrix<double>& C_alpha, const TransformationMatrix<double>& C_beta, const TransformationMatrix<double>& C, const QCMatrix<double>& S) const {


    // Make a reference copy in order to improve readibility of the following code.
    const auto& onv_of = *this;


    // Calculate the transformation matrices between both sets of spin-orbitals.
    TransformationMatrix<double> T_alpha = C.adjoint() * S * C_alpha;
    TransformationMatrix<double> T_beta = C.adjoint() * S * C_beta;


    // T's columns should be the ones occupied in the 'of'-ONV.
    // T's rows should be the ones occupied in the 'on'-ONV.
    // While waiting for Eigen 3.4 to release (which has better slicing APIs), we'll remove the UNoccupied rows/columns.
    const auto unoccupied_indices_of_alpha = onv_of.onv(Spin::alpha).unoccupiedIndices();
    const auto unoccupied_indices_on_alpha = onv_on.onv(Spin::alpha).unoccupiedIndices();

    T_alpha.removeColumns(unoccupied_indices_of_alpha);
    T_alpha.removeRows(unoccupied_indices_on_alpha);


    const auto unoccupied_indices_of_beta = onv_of.onv(Spin::beta).unoccupiedIndices();
    const auto unoccupied_indices_on_beta = onv_on.onv(Spin::beta).unoccupiedIndices();

    T_beta.removeColumns(unoccupied_indices_of_beta);
    T_beta.removeRows(unoccupied_indices_on_beta);


    // The requested overlap is the determinant of the product of the resulting smaller matrices.
    return (T_alpha * T_beta).determinant();
}


/**
 *  @param sigma                alpha or beta
 * 
 *  @return the number of sigma-electrons this spin-resolved ONV describes
 */
size_t SpinResolvedONV::numberOfElectrons(const Spin sigma) const {

    switch (sigma) {
    case Spin::alpha: {
        return this->onv_alpha.numberOfElectrons();
        break;
    }

    case Spin::beta: {
        return this->onv_beta.numberOfElectrons();
        break;
    }
    }
}


/**
 *  @param sigma                alpha or beta
 * 
 *  @return the number of sigma-spatial orbitals/spin-orbitals that this ONV is expressed with
 */
size_t SpinResolvedONV::numberOfSpatialOrbitals(const Spin sigma) const {

    switch (sigma) {
    case Spin::alpha: {
        return this->onv_alpha.numberOfSpinors();
        break;
    }

    case Spin::beta: {
        return this->onv_beta.numberOfSpinors();
        break;
    }
    }
}


/**
 *  @param sigma                alpha or beta
 * 
 *  @return the ONV that describes the occupations of the sigma-spin orbitals.
 */
const SpinUnresolvedONV& SpinResolvedONV::onv(const Spin sigma) const {

    switch (sigma) {

    case Spin::alpha: {
        return this->onv_alpha;
        break;
    }

    case Spin::beta: {
        return this->onv_beta;
        break;
    }
    }
}


}  // namespace GQCP