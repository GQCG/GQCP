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
 *  Calculate the overlap between this and another spin-resolved ONV expressed in R/U-spinor bases.
 * 
 *  @param other                        the other spin-resolved ONV
 *  @param this_spinor_basis            the restricted spin-orbital basis in which this ONV is expressed
 *  @param other_spinor_basis           the unrestricted spin-orbital basis in which the other ONV is expressed
 * 
 *  @param the overlap between this and the other spin-resolved ONV
 */
double SpinResolvedONV::calculateOverlap(const SpinResolvedONV& other, const RSpinorBasis<double, GTOShell>& this_spinor_basis, const USpinorBasis<double, GTOShell>& other_spinor_basis) const {

    const auto S = this_spinor_basis.overlap().parameters();  // the overlap matrix in AO basis
    const auto S_alpha = other_spinor_basis.overlap(Spin::alpha).parameters();
    const auto S_beta = other_spinor_basis.overlap(Spin::beta).parameters();
    if (!(S.isApprox(S_alpha, 1.0e-08)) || !(S.isApprox(S_beta, 1.0e-08))) {
        throw std::invalid_argument("SpinResolvedONV::calculateOverlap(const SpinResolvedONV&, const RSpinorBasis<double, GTOShell>&, const USpinorBasis<double, GTOShell>&): The given spinor bases are not expressed using the same scalar orbital basis.");
    }


    // Prepare some parameters.
    const auto& C = this_spinor_basis.coefficientMatrix();
    const auto& C_alpha = other_spinor_basis.coefficientMatrix(Spin::alpha);
    const auto& C_beta = other_spinor_basis.coefficientMatrix(Spin::beta);


    // Calculate the overlap matrices between the restricted and unrestricted spinor bases.
    MatrixX<double> T_alpha = C.adjoint() * S * C_alpha;
    MatrixX<double> T_beta = C.adjoint() * S * C_beta;


    // Create the smaller 'occupied' matrices by deleting columns and rows belonging to unoccupied spin-orbitals.
    const auto unoccupied_indices_alpha = this->onv(Spin::alpha).unoccupiedIndices();
    T_alpha.removeRows(unoccupied_indices_alpha);
    T_alpha.removeColumns(unoccupied_indices_alpha);

    const auto unoccupied_indices_beta = this->onv(Spin::beta).unoccupiedIndices();
    T_beta.removeRows(unoccupied_indices_beta);
    T_beta.removeColumns(unoccupied_indices_beta);


    // The overlap is the determinant of the product of the resulting matrices.
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