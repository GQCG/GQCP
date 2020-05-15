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
#include "Mathematical/Representation/BlockMatrix.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


namespace GQCP {


/**
 *  The coupled-cluster T1 amplitudes t_i^a. According to context, this class may either represent restricted (i.e. spatial-orbital) amplitudes, or generalized (spinor) amplitudes.
 */
template <typename _Scalar>
class T1Amplitudes {
public:
    using Scalar = _Scalar;


private:
    OrbitalSpace orbital_space;  // the orbital space which covers the occupied-virtual separation

    BlockMatrix<Scalar> t;  // the T1-amplitudes as a block matrix, implementing easy operator(i,a) calls


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Construct T1-amplitudes given their representation as a BlockMatrix and and explicit occupied-virtual orbital space.
     * 
     *  @param t                            the T1-amplitudes as a block matrix, implementing easy operator(i,a) calls
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation, indicating which indices in the block matrix are occupied and which are virtual
     */
    T1Amplitudes(const BlockMatrix<double>& t, const OrbitalSpace& orbital_space) :
        orbital_space {orbital_space},
        t {t} {}


    /**
     *  Construct the T1-amplitudes given their representation as a BlockMatrix and an implicit occupied-virtual orbital space determined by the given number of occupied orbitals and total number of orbitals.
     * 
     *  @param t                the T1-amplitudes as a block matrix, implementing easy operator(i,a) calls
     *  @param N                the number of occupied orbitals
     *  @param M                the total number of orbitals
     */
    T1Amplitudes(const BlockMatrix<double>& t, const size N, const size_t M) :
        T1Amplitudes(t, OrbitalSpace::OccupiedVirtual(N, M)) {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Create perturbative T1-amplitudes using an explicit orbital space.
     * 
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal spinor basis
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation
     * 
     *  @return T1-amplitudes calculated from an initial perturbative result
     */
    static T1Amplitudes<Scalar> Perturbative(const SQHamiltonian<Scalar>& sq_hamiltonian, const OrbitalSpace& orbital_space) {

        // Zero-initialize a BlockMatrix of the appropriate dimensions.
        const auto N = orbital_space.numberOfOccupiedOrbitals();
        const auto M = orbital_space.numberOfVirtualOrbitals();
        BlockMatrix<double> t1 {0, N, N, M};  // a block matrix suitable for occupied-virtual objects, like these T1-amplitudes t_i^a


        // Provide the perturbative T1-amplitudes.
        const auto F = sq_hamiltonian.calculateInactiveFockian(orbital_space);

        for (const auto& i : orbital_space.occupiedIndices()) {
            for (const auto& a : orbital_space.virtualIndices()) {
                const auto denominator = F(i, i) - F(a, a);

                t1(i, a) = F(a, i) / denominator;
            }
        }

        return T1Amplitudes<Scalar>(t1, orbital_space);
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
     *  @param a            a virtual index
     * 
     *  @return the T1-amplitude corresponding to t_i^a
     */
    double operator()(const size_t i, const size_t a) const { return this->t(i, a); }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the T1-amplitudes as a BlockMatrix
     */
    const BlockMatrix<double>& asBlockMatrix() const { return this->t; }

    /**
     *  @return the orbital space for these T1-amplitudes, which covers the occupied-virtual separation
     */
    const OrbitalSpace& orbitalSpace() const { return this->orbital_space; }
};


}  // namespace GQCP
