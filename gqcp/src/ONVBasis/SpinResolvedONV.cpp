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
 *  Create a spin-resolved ONV that represents the RHF single Slater determinant.
 * 
 *  @param K            the number of spatial orbitals
 *  @param N_P          the number of electron pairs
 * 
 *  @param a spin-resolved ONV that represents the RHF single Slater determinant
 */
SpinResolvedONV SpinResolvedONV::RHF(const size_t K, const size_t N_P) {

    const auto sigma_onv = SpinUnresolvedONV::GHF(K, N_P);  // for sigma = alpha and sigma = beta

    return SpinResolvedONV(sigma_onv, sigma_onv);
}


/**
 *  Create a spin-resolved ONV that represents the UHF single Slater determinant.
 * 
 *  @param K            the number of spatial orbitals
 *  @param N_P          the number of electron pairs
 * 
 *  @param a spin-resolved ONV that represents the UHF single Slater determinant
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
 *  Calculate the overlap <on|of>: the projection of between this spin-resolved ONV ('of') and another spin-resolved ONV ('on'), expressed in different R/U-spinor bases.
 * 
 *  @param onv_on                       the spin-resolved ONV that should be projected on
 *  @param spinor_basis_of              the unrestricted spin-orbital basis in which this ONV (the 'of'-ONV) is expressed
 *  @param spinor_basis_on              the restricted spin-orbital basis in which the 'on'-ONV is expressed.
 * 
 *  @return the overlap element <on|of>
 * 
 *  @note This method can be used to project UHF-ONVs onto RHF-ONVs, by calling
 *          uhf_onv.calculateProjection(rhf_onv, USpinorBasis, RSpinorBasis)
 */
double SpinResolvedONV::calculateProjection(const SpinResolvedONV& onv_on, const USpinorBasis<double, GTOShell>& spinor_basis_of, const RSpinorBasis<double, GTOShell>& spinor_basis_on) const {


    // Determine the overlap matrices of the underlying scalar orbital bases, which is needed later on.
    auto S = spinor_basis_on.overlap().parameters();                         // the overlap matrix of the restricted MOs/spin-orbitals
    S.basisTransformInPlace(spinor_basis_on.coefficientMatrix().inverse());  // now in AO basis

    auto S_alpha = spinor_basis_of.overlap(Spin::alpha).parameters();                         // the overlap matrix of the alpha spin-orbitals
    S_alpha.basisTransformInPlace(spinor_basis_of.coefficientMatrix(Spin::alpha).inverse());  // now in AO basis

    auto S_beta = spinor_basis_of.overlap(Spin::beta).parameters();                         // the overlap matrix of the beta spin-orbitals
    S_beta.basisTransformInPlace(spinor_basis_of.coefficientMatrix(Spin::beta).inverse());  // now in AO basis

    if (!(S.isApprox(S_alpha, 1.0e-08)) || !(S.isApprox(S_beta, 1.0e-08))) {
        throw std::invalid_argument("SpinResolvedONV::calculateOverlap(const SpinResolvedONV&, const RSpinorBasis<double, GTOShell>&, const USpinorBasis<double, GTOShell>&): The given spinor bases are not expressed using the same scalar orbital basis.");
    }


    // Prepare some parameters.
    const auto onv_of = *this;
    const auto& C = spinor_basis_on.coefficientMatrix();
    const auto& C_alpha = spinor_basis_of.coefficientMatrix(Spin::alpha);
    const auto& C_beta = spinor_basis_of.coefficientMatrix(Spin::beta);


    // Calculate the transformation matrices between both sets of spin-orbitals.
    TransformationMatrix<double> T_alpha = C.adjoint() * S * C_alpha;
    TransformationMatrix<double> T_beta = C.adjoint() * S * C_beta;


    // T's columns should be the ones occupied in the 'of'-ONV.
    // T's rows should be the ones occupied in the 'on'-ONV.
    // While waiting for Eigen 3.4. to release (which has better slicing APIs), we'll remove the UNoccupied rows/columns.
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
size_t SpinResolvedONV::numberOfElectrons(Spin sigma) const {

    switch (sigma) {
    case Spin::alpha: {
        return this->onv_alpha.numberOfElectrons();
    }

    case Spin::beta: {
        return this->onv_beta.numberOfElectrons();
    }
    }
}


/**
 *  @param sigma                alpha or beta
 * 
 *  @return the number of sigma-spatial orbitals/spin-orbitals that this ONV is expressed with
 */
size_t SpinResolvedONV::numberOfSpatialOrbitals(Spin sigma) const {

    switch (sigma) {
    case Spin::alpha: {
        return this->onv_alpha.numberOfSpinors();
    }

    case Spin::beta: {
        return this->onv_beta.numberOfSpinors();
    }
    }
}


/**
 *  @param sigma                alpha or beta
 * 
 *  @return the ONV that describes the occupations of the sigma-spin orbitals.
 */
const SpinUnresolvedONV& SpinResolvedONV::onv(Spin sigma) const {

    switch (sigma) {

    case Spin::alpha: {
        return this->onv_alpha;
    }

    case Spin::beta: {
        return this->onv_beta;
    }
    }
}


}  // namespace GQCP