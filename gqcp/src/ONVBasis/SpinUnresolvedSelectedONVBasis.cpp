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

#include "ONVBasis/SpinUnresolvedSelectedONVBasis.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  Construct an empty spin-unresolved selected ONV basis.
 *
 *  @param M            The number of spinors.
 *  @param N            The number of electrons.
 */
SpinUnresolvedSelectedONVBasis::SpinUnresolvedSelectedONVBasis(const size_t M, const size_t N) :
    M {M},
    N {N} {}


/**
 *  Generate a `SpinUnresolvedSelectedONVBasis` from a full spin-unresolved ONV basis.
 *
 *  @param onv_basis        The full spin-unresolved ONV basis.
 */
SpinUnresolvedSelectedONVBasis::SpinUnresolvedSelectedONVBasis(const SpinUnresolvedONVBasis& onv_basis) :
    SpinUnresolvedSelectedONVBasis(onv_basis.numberOfOrbitals(), onv_basis.numberOfElectrons()) {

    // Loop through the full spin-unresolved ONV basis and insert every ONV.
    std::vector<SpinUnresolvedONV> onvs;

    onv_basis.forEach([&onvs](const SpinUnresolvedONV& onv, const size_t I) {
        onvs.push_back(onv);
    });

    this->onvs = onvs;
}


/*
 *  MARK: Modifying
 */

/**
 *  Expand this ONV basis with the given spin-unresolved ONV.
 * 
 *  @param onv          The ONV that should be included in this ONV basis.
 */
void SpinUnresolvedSelectedONVBasis::expandWith(const SpinUnresolvedONV& onv) {

    if (onv.numberOfElectrons() != this->numberOfElectrons()) {
        throw std::invalid_argument("SpinUnresolvedSelectedONVBasis::expandWith(const SpinUnesolvedONV&): The given ONV's number of electrons is not compatible with the number of electrons for this ONV basis.");
    }

    if (onv.numberOfSpinors() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinUnresolvedSelectedONVBasis::expandWith(const std::vector<SpinResolvedONV>&): The given ONV's number of orbitals is not compatible with the number of orbitals for this ONV basis.");
    }

    this->onvs.push_back(onv);
}


/**
 *  Expand this ONV basis with the given spin-unresolved ONVs.
 * 
 *  @param onvs         The ONVs that should be included in this ONV basis.
 */
void SpinUnresolvedSelectedONVBasis::expandWith(const std::vector<SpinUnresolvedONV>& onvs) {

    for (const auto& onv : onvs) {
        this->onvs.push_back(onv);
    }
}


}  // namespace GQCP
