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


#include "DensityMatrix/SpinResolved1DM.hpp"
#include "DensityMatrix/SpinResolved2DM.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


namespace GQCP {


/**
 *  Calculate the expectation value of the square of the spin angular momentum operator.
 * 
 *  @tparam Scalar        The scalar type of the density matrices.
 * 
 *  @param one_DMs        The one-electron density matrix.
 *  @param two_DMs        The two-electron density matrix.
 *
 *  @return The expectation value of the square of the spin angular momentum operator.
 */
template <typename Scalar>
double calculateSpinSquared(const SpinResolved1DM<Scalar>& one_DMs, const SpinResolved2DM<Scalar>& two_DMs) {

    double sz = calculateSpinZ(one_DMs);
    double s_squared = -sz;
    const size_t K = one_DMs.alpha().numberOfOrbitals();
    for (size_t p = 0; p < K; p++) {
        s_squared += one_DMs.alpha().matrix()(p, p);                                        // One-electron partition of S+S_
        s_squared += (one_DMs.alpha().matrix()(p, p) + one_DMs.beta().matrix()(p, p)) / 4;  // One-electron partition of S^2
        for (size_t q = 0; q < K; q++) {
            s_squared += -two_DMs.alphaBeta().tensor()(p, q, q, p);                                                                                                                                        // Two-electron partition  S+S_ .
            s_squared += (two_DMs.alphaAlpha().tensor()(p, p, q, q) + two_DMs.betaBeta().tensor()(p, p, q, q) - two_DMs.alphaBeta().tensor()(p, p, q, q) - two_DMs.betaAlpha().tensor()(p, p, q, q)) / 4;  // Two-electron partition of S^2 .
        }
    }
    return s_squared;
}


/**
 *  Calculate the expectation value of the z-component of the spin angular momentum operator.
 * 
 *  @tparam Scalar         The scalar type of the density matrix.
 * 
 *  @param one_DMs         The one-electron density matrix.
 *
 *  @return The expectation value of spin in the z-direction.
 */
template <typename Scalar>
double calculateSpinZ(const SpinResolved1DM<Scalar>& one_DMs) {
    return one_DMs.spinDensity().matrix().trace() / 2;
}


}  // namespace GQCP
