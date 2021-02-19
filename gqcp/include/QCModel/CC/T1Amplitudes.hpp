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
#include "QCModel/CC/SinglesAmplitudes.hpp"


namespace GQCP {


/**
 *  The coupled-cluster T1 amplitudes t_i^a. According to context, this class may either represent restricted (i.e. spatial-orbital) amplitudes, or generalized (spinor) amplitudes.
 */
template <typename _Scalar>
class T1Amplitudes:
    public SinglesAmplitudes<_Scalar, T1Amplitudes<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Self = T1Amplitudes<Scalar>;


public:
    /*
     *  CONSTRUCTORS
     */

    // Inherit `SinglesAmplitudes`' constructors.
    using SinglesAmplitudes<Scalar, T1Amplitudes<Scalar>>::SinglesAmplitudes;


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Create perturbative T1-amplitudes using an explicit orbital space.
     * 
     *  @param f                            the (inactive) Fock matrix
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation
     * 
     *  @return T1-amplitudes calculated from an initial perturbative result
     */
    static T1Amplitudes<Scalar> Perturbative(const SquareMatrix<Scalar>& f, const OrbitalSpace& orbital_space) {

        // Zero-initialize a matrix representation for the (occupied-virtual) T1-amplitudes t_i^a.
        auto t1 = orbital_space.initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual);

        // Provide the perturbative T1-amplitudes, by setting the RHS amplitudes in Stanton1991, equation (1) to zero.
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                const auto denominator = f(i, i) - f(a, a);

                t1(i, a) = f(i, a) / denominator;
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
    static T1Amplitudes<Scalar> Perturbative(const RSQHamiltonian<Scalar>& sq_hamiltonian, const size_t N, const size_t M) {

        // Create the implicit orbital space for N occupied orbitals and M total orbitals.
        const auto orbital_space = OrbitalSpace::Implicit({{OccupationType::k_occupied, N}, {OccupationType::k_virtual, M - N}});

        return T1Amplitudes<Scalar>::Perturbative(sq_hamiltonian, orbital_space);
    }
};

}  // namespace GQCP
