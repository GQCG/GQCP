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


#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "Mathematical/Representation/BlockRankFourTensor.hpp"


namespace GQCP {


/**
 *  The coupled-cluster T2-amplitudes t_{ij}^{ab}. According to context, this class may either represent restricted (i.e. spatial-orbital) amplitudes, or generalized (spinor) amplitudes.
 */
template <typename _Scalar>
class T2Amplitudes {
public:
    using Scalar = _Scalar;

private:
    OrbitalSpace orbital_space;  // the orbital space which covers the occupied-virtual separation

    BlockRankFourTensor<Scalar> t;  // the T2-amplitudes as a block rank-four tensor, implementing easy operator(i,j,a,b) calls


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Construct T2-amplitudes given their representation as a BlockRankFourTensor and and explicit occupied-virtual orbital space.
     * 
     *  @param t                            the T2-amplitudes as a block matrix, implementing easy operator(i,j,a,b) calls
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation, indicating which indices in the block matrix are occupied and which are virtual
     */
    T1Amplitudes(const BlockRankFourTensor<double>& t, const OrbitalSpace& orbital_space) :
        orbital_space {orbital_space},
        t {t} {}


    /**
     *  Construct the T2-amplitudes given their representation as a BlockRankFourTensor and an implicit occupied-virtual orbital space determined by the given number of occupied orbitals and total number of orbitals.
     * 
     *  @param t                the T2-amplitudes as a block matrix, implementing easy operator(i,j,a,b) calls
     *  @param N                the number of occupied orbitals
     *  @param M                the total number of orbitals
     */
    T1Amplitudes(const BlockMatrix<double>& t, const size N, const size_t M) :
        T2Amplitudes(t, OrbitalSpace::OccupiedVirtual(N, M)) {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Create perturbative T2-amplitudes using an explicit orbital space.
     * 
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal spinor basis
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation
     * 
     *  @return T2-amplitudes calculated from an initial perturbative result
     */
    static T2Amplitudes<Scalar> Perturbative(const SQHamiltonian<Scalar>& sq_hamiltonian, const OrbitalSpace& orbital_space) {

        // Zero-initialize a BlockMatrix of the appropriate dimensions.
        const auto N = orbital_space.numberOfOccupiedOrbitals();
        const auto M = orbital_space.numberOfVirtualOrbitals();
        BlockRankFourTensor<double> t2 {0, N, 0, N,
                                        N, M, N, M};  // a block rank-four tensor suitable for occupied-occupied-virtual-virtual objects, like these T2-amplitudes t_{ij}^{ab}

        // Provide the perturbative T2-amplitudes.
        const auto F = sq_hamiltonian.calculateInactiveFockian(orbital_space);

        const auto& g_chemists = sq_hamiltonian.twoElectron().parameters();
        const auto V_A = g_chemists.convertedToPhysicistsNotation().antisymmetrized();

        for (const auto& i : orbital_space.occupiedIndices()) {
            for (const auto& j : orbital_space.occupiedIndices()) {
                for (const auto& a : orbital_space.virtualIndices()) {
                    for (const auto& b : orbital_space.virtualIndices()) {
                        const auto denominator = F(i, i) + F(j, j) - F(a, a) - F(b, b);

                        t2(i, j, a, b) = V_A(a, b, i, j) / denominator;
                    }
                }
            }
        }

        return T2Amplitudes<Scalar>(t2, orbital_space);
    }


    /**
     *  Create perturbative T1-amplitudes using an implicit occupied-virtual orbital space determined by the given number of occupied orbitals and total number of orbitals.
     * 
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal spinor basis
     *  @param N                            the number of occupied orbitals
     *  @param M                            the total number of orbitals
     * 
     *  @return T1-amplitudes calculated from an initial perturbative result
     */
    static T1Amplitudes<Scalar> Perturbative(const SQHamiltonian<Scalar>& sq_hamiltonian, const size_t N, const size_t M) {

        // Create the implicit orbital space for N occupied orbitals and M total orbitals.
        const auto orbital_space = OrbitalSpace::OccupiedVirtual(N, M);

        return T1Amplitudes<Scalar>::Perturbative(sq_hamiltonian, orbital_space);
    }


    /*
     *  OPERATORS
     */

    /**
     *  @param i            an occupied index
     *  @param j            an occupied index
     *  @param a            a virtual index
     *  @param b            a virtual index
     * 
     *  @return the T2-amplitude corresponding to t_{ij}^{ab}
     */
    double operator()(const size_t i, const size_t j, const size_t a, const size_t b) const { return this->t(i, j, a, b); }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the T2-amplitudes as a BlockRankFourTensor
     */
    const BlockRankFourTensor<double>& asBlockRankFourTensor() const { return this->t; }


    /**
     *  @return the orbital space for these T2-amplitudes, which covers the occupied-virtual separation
     */
    const OrbitalSpace& orbitalSpace() const { return this->orbital_space; }
};


}  // namespace GQCP
