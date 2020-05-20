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

    SQHamiltonian<Scalar> sq_hamiltonian;  // the Hamiltonian expressed in an orthonormal spinor basis


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Initialize a CCSD algorithmic environment with given T1- and T2-amplitudes.
     * 
     *  @param t1_amplitudes            the initial T1-amplitudes
     *  @param t2_amplitudes            the initial T2-amplitudes
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal spinor basis
     */
    CCSDEnvironment(const T1Amplitudes<Scalar>& t1_amplitudes, const T2Amplitudes<Scalar>& t2_amplitudes, const SQHamiltonian<Scalar>& sq_hamiltonian) :
        electronic_energies {QCModel::CCSD<Scalar>::calculateCorrelationEnergy(sq_hamiltonian, t1_amplitudes, t2_amplitudes)},  // already calculate the initial CCSD energy correction
        t1_amplitudes {t1_amplitudes},
        t2_amplitudes {t2_amplitudes},
        sq_hamiltonian {sq_hamiltonian} {}


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
    static CCSDEnvironment<Scalar> Perturbative(const SQHamiltonian<Scalar>& sq_hamiltonian, const OrbitalSpace& orbital_space) {

        const auto t1_amplitudes = T1Amplitudes<Scalar>::Perturbative(sq_hamiltonian, orbital_space);
        const auto t2_amplitudes = T2Amplitudes<Scalar>::Perturbative(sq_hamiltonian, orbital_space);

        return CCSDEnvironment<Scalar>(t1_amplitudes, t2_amplitudes, sq_hamiltonian);
    }
};


}  // namespace GQCP
