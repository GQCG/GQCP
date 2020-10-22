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
#include "QCModel/CC/CCD.hpp"
#include "QCModel/CC/CCSD.hpp"
#include "QCModel/CC/T1Amplitudes.hpp"
#include "QCModel/CC/T2Amplitudes.hpp"

#include <deque>


namespace GQCP {


/**
 *  An algorithmic environment suitable for coupled-cluster calculations up to the CCSD level.
 * 
 *  @tparam _Scalar             the scalar type the amplitudes
 */
template <typename _Scalar>
class CCSDEnvironment {
public:
    using Scalar = _Scalar;


public:
    std::deque<double> electronic_energies;  // the electronic correlation energy

    std::deque<T1Amplitudes<Scalar>> t1_amplitudes;
    std::deque<T2Amplitudes<Scalar>> t2_amplitudes;

    SquareMatrix<Scalar> f;        // the (inactive) Fock matrix
    QCRankFourTensor<Scalar> V_A;  // the antisymmetrized two-electron integrals (in physicist's notation)

    ImplicitMatrixSlice<Scalar> F1;  // represents equation (3) in Stanton1991
    ImplicitMatrixSlice<Scalar> F2;  // represents equation (4) in Stanton1991
    ImplicitMatrixSlice<Scalar> F3;  // F3 represents equation (5) in Stanton1991

    ImplicitRankFourTensorSlice<Scalar> W1;  // represents equation (6) in Stanton1991
    ImplicitRankFourTensorSlice<Scalar> W2;  // represents equation (7) in Stanton1991
    ImplicitRankFourTensorSlice<Scalar> W3;  // represents equation (8) in Stanton1991

    ImplicitRankFourTensorSlice<Scalar> tau2;        // represents equation (10) in Stanton1991
    ImplicitRankFourTensorSlice<Scalar> tau2_tilde;  // represents equation (9) in Stanton1991


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Initialize a CCSD algorithmic environment with given T1- and T2-amplitudes.
     * 
     *  @param t1_amplitudes            the initial T1-amplitudes
     *  @param t2_amplitudes            the initial T2-amplitudes
     *  @param f                        the (inactive) Fock matrix
     *  @param V_A                      the antisymmetrized two-electron integrals (in physicist's notation)
     */
    CCSDEnvironment(const T1Amplitudes<Scalar>& t1_amplitudes, const T2Amplitudes<Scalar>& t2_amplitudes, const SquareMatrix<Scalar>& f, const QCRankFourTensor<Scalar>& V_A) :
        electronic_energies {QCModel::CCSD<Scalar>::calculateCorrelationEnergy(f, V_A, t1_amplitudes, t2_amplitudes)},  // already calculate the initial CCSD energy correction
        t1_amplitudes {t1_amplitudes},
        t2_amplitudes {t2_amplitudes},
        f {f},
        V_A {V_A} {}

    /**
     *  Initialize a CCD algorithmic environment with given T2-amplitudes.
     * 
     *  @param t2_amplitudes            the initial T2-amplitudes
     *  @param f                        the (inactive) Fock matrix
     *  @param V_A                      the antisymmetrized two-electron integrals (in physicist's notation)
     */
    CCSDEnvironment(const T2Amplitudes<Scalar>& t2_amplitudes, const SquareMatrix<Scalar>& f, const QCRankFourTensor<Scalar>& V_A) :
        electronic_energies {QCModel::CCD<Scalar>::calculateCorrelationEnergy(f, V_A, t2_amplitudes)},  // already calculate the initial CCD energy correction
        t2_amplitudes {t2_amplitudes},
        f {f},
        V_A {V_A} {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Initialize a CCSD algorithmic environment with initial guesses for the T1- and T2-amplitudes based on perturbation theory.
     * 
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal spinor basis
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation
     * 
     *  @return an algorithmic environment suitable for coupled-cluster calculations up to the CCSD level.
     */
    static CCSDEnvironment<Scalar> PerturbativeCCSD(const GSQHamiltonian<Scalar>& sq_hamiltonian, const OrbitalSpace& orbital_space) {

        // For the CCSD environment equation, we need the inactive Fock matrix and the anti-symmetrized two-electron integrals in physicist's notation.
        const auto f = sq_hamiltonian.calculateInactiveFockian(orbital_space).parameters();

        const auto& g_chemists = sq_hamiltonian.twoElectron().parameters();
        const auto V_A = g_chemists.convertedToPhysicistsNotation().antisymmetrized();


        const auto t1_amplitudes = T1Amplitudes<Scalar>::Perturbative(f, orbital_space);
        const auto t2_amplitudes = T2Amplitudes<Scalar>::Perturbative(f, V_A, orbital_space);

        return CCSDEnvironment<Scalar>(t1_amplitudes, t2_amplitudes, f, V_A);
    }

    /**
     *  Initialize a CCD algorithmic environment with initial guesses for the T2-amplitudes based on perturbation theory.
     * 
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal spinor basis
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation
     * 
     *  @return an algorithmic environment suitable for CCD calculations.
     */
    static CCSDEnvironment<Scalar> PerturbativeCCD(const GSQHamiltonian<Scalar>& sq_hamiltonian, const OrbitalSpace& orbital_space) {

        // For the CCSD environment equation, we need the inactive Fock matrix and the anti-symmetrized two-electron integrals in physicist's notation.
        const auto f = sq_hamiltonian.calculateInactiveFockian(orbital_space).parameters();

        const auto& g_chemists = sq_hamiltonian.twoElectron().parameters();
        const auto V_A = g_chemists.convertedToPhysicistsNotation().antisymmetrized();

        const auto t2_amplitudes = T2Amplitudes<Scalar>::Perturbative(f, V_A, orbital_space);

        return CCSDEnvironment<Scalar>(t2_amplitudes, f, V_A);
    }
};


}  // namespace GQCP
