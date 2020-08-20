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


#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/DensityMatrices/SpinResolvedOneDM.hpp"
#include "Processing/DensityMatrices/SpinResolvedTwoDM.hpp"


namespace GQCP {


/**
 *  Calculate the expectation value of the square of the spin angular momentum operator
 * 
 *  @tparam Scalar         the scalar type of the density matrices
 * 
 *  @param one_rdms        all the one-electron density matrices
 *  @param two_rdms        all the two-electron density matrices
 *
 *  @return the expectation value of the square of the spin angular momentum operator
 */
template <typename Scalar>
double calculateSpinSquared(const SpinResolvedOneDM<Scalar>& one_rdms, const SpinResolvedTwoDM<Scalar>& two_rdms) {

    double sz = calculateSpinZ(one_rdms);
    double s_squared = -sz;
    const size_t K = one_rdms.dimension();
    for (size_t p = 0; p < K; p++) {
        s_squared += one_rdms.one_rdm_aa(p, p);                                    // One-electron partition of S+S_
        s_squared += (one_rdms.one_rdm_aa(p, p) + one_rdms.one_rdm_bb(p, p)) / 4;  // One-electron partition of S^2
        for (size_t q = 0; q < K; q++) {
            s_squared += -two_rdms.two_rdm_aabb(p, q, q, p);                                                                                                                   // Two-electron partition  S+S_
            s_squared += (two_rdms.two_rdm_aaaa(p, p, q, q) + two_rdms.two_rdm_bbbb(p, p, q, q) - two_rdms.two_rdm_aabb(p, p, q, q) - two_rdms.two_rdm_bbaa(p, p, q, q)) / 4;  // Two-electron partition of S^2
        }
    }
    return s_squared;
}


/**
 *  Calculate the expectation value of the z-component of the spin angular momentum operator
 * 
 *  @tparam Scalar          the scalar type
 * 
 *  @param one_rdms         all the one-electron density matrices
 *
 *  @return expectation value of spin in the z direction
 */
template <typename Scalar>
double calculateSpinZ(const SpinResolvedOneDM<Scalar>& one_rdms) {
    return one_rdms.spinDensityRDM().trace() / 2;
}


}  // namespace GQCP
