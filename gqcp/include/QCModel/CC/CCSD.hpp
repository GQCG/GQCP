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
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/CC/T1Amplitudes.hpp"
#include "QCModel/CC/T2Amplitudes.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The CCSD (coupled-cluster singles and doubles) wave function model.
 * 
 *  @tparam _Scalar             the scalar type of the amplitudes
 */
template <typename _Scalar>
class CCSD {
public:
    using _Scalar = Scalar;


private:
public:
    /*
     *  STATIC PUBLIC METHODS
     */

    /**
     *  Calculate the CCSD correlation energy.
     * 
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal spinor basis
     *  @param t1                           the T1-amplitudes
     *  @param t2                           the T2-amplitudes
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation
     * 
     *  @return the CCSD correlation energy
     */
    static double calculateCorrelationEnergy(const SQHamiltonian<Scalar>& sq_hamiltonian, const T1Amplitudes<Scalar>& t1, const T2Amplitudes<Scalar>& t2, const OrbitalSpace& orbital_space) {

        // For the CCSD energy equation, we need the inactive Fock matrix and the anti-symmetrized two-electron integrals in physicist's notation.
        const auto F = sq_hamiltonian.calculateInactiveFockian(orbital_space);

        const auto& g_chemists = sq_hamiltonian.twoElectron().parameters();
        const auto V_A = g_chemists.convertedToPhysicistsNotation().antisymmetrized();


        // A KISS implementation of the CCSD energy correction.
        // The implementation is in line with Crawford2000 "Chapter 2: An Introduction to Coupled Cluster Theory for Computational Chemists", eq. [134].
        double E = 0.0;
        for (const auto& i : orbital_space.occupiedIndices()) {
            for (const auto& a : orbital_space.virtualIndices()) {
                E += F(i, a) * t1(i, a);
            }
        }

        for (const auto& i : orbital_space.occupiedIndices()) {
            for (const auto& j : orbital_space.occupiedIndices()) {
                for (const auto& a : orbital_space.virtualIndices()) {
                    for (const auto& b : orbital_space.virtualIndices()) {
                        E += 0.25 * V_A(i, j, a, b) * t2(i, j, a, b);

                        E += 0.5 * V_A(i, j, a, b) * t1(i, a) * t1(j, b);
                    }
                }
            }
        }

        return E;
    }


    /**
     *  Calculate the value for one of the CCSD T1-amplitude equations, evaluated at the given T1- and T2-amplitudes
     *      f_i^a = <Phi_i^a| H |Phi_0>             with H the similarity-transformed normal-ordered Hamiltonian
     * 
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal spinor basis
     *  @param t1                           the T1-amplitudes
     *  @param t2                           the T2-amplitudes
     *  @param i                            the (occupied) subscript for the amplitude equation
     *  @param a                            the (virtual) superscript for the amplitude equation
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation
     * 
     *  @return the value for one of the CCSD T1-amplitude equations
     */
    static double calculateT1AmplitudeEquation(const SQHamiltonian<double>& sq_hamiltonian, const T1Amplitudes<Scalar>& t1, const T2Amplitudes<Scalar>& t2, const size_t i, const size_t a, const OrbitalSpace& orbital_space) {

        // For the CCSD T1-amplitude equations equations, we need the inactive Fock matrix and the anti-symmetrized two-electron integrals in physicist's notation.
        const auto F = sq_hamiltonian.calculateInactiveFockian(orbital_space).parameters();

        const auto& g_chemists = sq_hamiltonian.twoElectron().parameters();
        const auto V_A = g_chemists.convertedToPhysicistsNotation().antisymmetrized();


        // A KISS implementation of the CCSD T1-amplitude equations.
        // The implementation is in line with Crawford2000 "Chapter 2: An Introduction to Coupled Cluster Theory for Computational Chemists", eq. [152].
        double value = 0.0;

        value += F(a, i);

        for (const auto& c : orbital_space.virtualIndices()) {
            value += F(a, c) * t1(i, c);
        }

        for (const auto& k : orbital_space.occupiedIndices()) {
            value -= F(k, i) * t1(k, a);
        }

        for (const auto& k : orbital_space.occupiedIndices()) {
            for (const auto& c : orbital_space.virtualIndices()) {
                value += V_A(k, a, c, i) * t1(k, c);

                value += F(k, c) * t2(i, k, a, c);

                value -= F(k, c) * t1(i, c) * t1(k, a);
            }
        }

        for (const auto& k : orbital_space.occupiedIndices()) {
            for (const auto& c : orbital_space.virtualIndices()) {
                for (const auto& d : orbital_space.virtualIndices()) {
                    value += 0.5 * V_A(k, a, c, d) * t2(k, i, c, d);

                    value -= V_A(k, a, c, d) * t1(k, c) * t1(i, d);
                }
            }
        }

        for (const auto& k : orbital_space.occupiedIndices()) {
            for (const auto& l : orbital_space.occupiedIndices()) {
                for (const auto& c : orbital_space.virtualIndices()) {
                    value -= 0.5 * V_A(k, l, c, i) * t2(k, l, c, a);

                    value -= V_A(k, l, c, i) * t1(k, c) * t1(l, a);
                }
            }
        }

        for (const auto& k : orbital_space.occupiedIndices()) {
            for (const auto& l : orbital_space.occupiedIndices()) {
                for (const auto& c : orbital_space.virtualIndices()) {
                    for (const auto& d : orbital_space.virtualIndices()) {
                        value -= V_A(k, l, c, d) * t1(k, c) * t1(i, d) * t1(l, a);

                        value += V_A(k, l, c, d) * t1(k, c) * t2(l, i, d, a);

                        value -= 0.5 * V_A(k, l, c, d) * t2(k, i, c, d) * t1(l, a);

                        value -= 0.5 * V_A(k, l, c, d) * t2(k, l, c, a) * t1(i, d);
                    }
                }
            }
        }


        return value;
    }


    /**
     *  Calculate the value for one of the CCSD T1-amplitude equations, evaluated at the given T1- and T2-amplitudes
     *      f_{ij}^{ab} = <Phi_{ij}^{ab}| H |Phi_0>             with H the similarity-transformed normal-ordered Hamiltonian
     * 
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal spinor basis
     *  @param t1                           the T1-amplitudes
     *  @param t2                           the T2-amplitudes
     *  @param i                            the (occupied) subscript for the amplitude equation
     *  @param j                            the other (occupied) subscript for the amplitude equation
     *  @param a                            the (virtual) superscript for the amplitude equation
     *  @param b                            the other (virtual) superscript for the amplitude equation
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation
     * 
     *  @return the value for one of the CCSD T1-amplitude equations
     */
    static double calculateT2AmplitudeEquation(const SQHamiltonian<double>& sq_hamiltonian, const T1Amplitudes<Scalar>& t1, const T2Amplitudes<Scalar>& t2, const size_t i, const size_t j, const size_t a, const size_t b, const OrbitalSpace& orbital_space) {


        // For the CCSD T2-amplitude equations equations, we need the inactive Fock matrix and the anti-symmetrized two-electron integrals in physicist's notation.
        const auto F = sq_hamiltonian.calculateInactiveFockian(orbital_space).parameters();

        const auto& g_chemists = sq_hamiltonian.twoElectron().parameters();
        const auto V_A = g_chemists.convertedToPhysicistsNotation().antisymmetrized();


        // A KISS implementation of the CCSD T2-amplitude equations.
        // The implementation is in line with Crawford2000 "Chapter 2: An Introduction to Coupled Cluster Theory for Computational Chemists", eq. [152].
        // In the following, P(pq) stands for the antisymmetric (i.e. with sign -1) permutation of the indices p and q
        double value = 0.0;

        value += V_A(a, b, i, j);

        for (const auto& c : orbital_space.virtualIndices()) {
            value += F(b, c) * t2(i, j, a, c);
            value -= F(a, c) * t2(i, j, b, c);

            value += V_A(a, b, c, j) * t1(i, c);
            value -= V_A(a, b, c, i) * t1(j, c);  // P(ij) applied
        }

        for (const auto& k : orbital_space.occupiedIndices()) {
            value -= F(k, j) * t2(i, k, a, b);
            value += F(k, i) * t2(j, k, a, b);

            value -= V_A(k, b, i, j) * t1(k, a);
            value += V_A(k, a, i, j) * t1(k, b);  // P(ab) applied
        }

        for (const auto& k : orbital_space.occupiedIndices()) {
            for (const auto& l : orbital_space.occupiedIndices()) {
                value += 0.5 * V_A(k, l, i, j) * t2(k, l, a, b);

                value += 0.5 * V_A(k, l, i, j) * t1(k, a) * t1(l, b);
                value -= 0.5 * V_A(k, l, i, j) * t1(k, b) * t1(l, a);  // P(ab) applied
            }
        }

        for (const auto& c : orbital_space.virtualIndices()) {
            for (const auto& d : orbital_space.virtualIndices()) {
                value += 0.5 * V_A(a, b, c, d) * t2(i, j, c, d);

                value += 0.5 * V_A(a, b, c, d) * t1(i, c) * t1(j, d);
                value -= 0.5 * V_A(a, b, c, d) * t1(j, c) * t1(i, d);  // P(ij) applied
            }
        }

        for (const auto& k : orbital_space.occupiedIndices()) {
            for (const auto& c : orbital_space.virtualIndices()) {
                value += V_A(k, b, c, j) * t2(i, k, a, c);
                value -= V_A(k, b, c, i) * t2(j, k, a, c);  // P(ij) applied
                value -= V_A(k, a, c, j) * t2(i, k, b, c);  // P(ab) applied
                value += V_A(k, a, c, i) * t2(j, k, b, c);  // P(ij) P(ab) applied

                value -= V_A(k, b, i, c) * t1(k, a) * t1(j, c);
                value += V_A(k, b, j, c) * t1(k, a) * t1(i, c);  // P(ij) applied
                value += V_A(k, a, i, c) * t1(k, b) * t1(j, c);  // P(ab) applied
                value -= V_A(k, a, j, c) * t1(k, b) * t1(i, c);  // P(ij) P(ab) applied

                value += F(k, c) * t1(k, a) * t2(i, j, b, c);
                value -= F(k, c) * t1(k, b) * t2(i, j, a, c);  // P(ab) applied

                value += F(k, c) * t1(i, c) * t2(j, k, a, b);
                value -= F(k, c) * t1(j, c) * t2(i, k, a, b);  // P(ij) applied
            }
        }

        for (const auto& k : orbital_space.occupiedIndices()) {
            for (const auto& l : orbital_space.occupiedIndices()) {
                for (const auto& c : orbital_space.virtualIndices()) {
                    value -= V_A(k, l, c, i) * t1(k, c) * t2(l, j, a, b);
                    value += V_A(k, l, c, j) * t1(k, c) * t2(l, i, a, b);  // P(ij) applied

                    value += V_A(k, l, i, c) * t1(l, a) * t2(j, k, b, c);
                    value -= V_A(k, l, j, c) * t1(l, a) * t2(i, k, b, c);  // P(ij) applied
                    value -= V_A(k, l, i, c) * t1(l, b) * t2(j, k, a, c);  // P(ab) applied
                    value += V_A(k, l, j, c) * t1(l, b) * t2(i, k, a, c);  // P(ij) P(ab) applied

                    value += 0.5 * V_A(k, l, c, j) * t1(i, c) * t2(k, l, a, b);
                    value -= 0.5 * V_A(k, l, c, i) * t1(j, c) * t2(k, l, a, b);  // P(ij) applied

                    value += 0.5 * V_A(k, l, c, j) * t1(i, c) * t1(k, a) * t1(l, b);
                    value -= 0.5 * V_A(k, l, c, i) * t1(j, c) * t1(k, a) * t1(l, b);  // P(ij) applied
                    value -= 0.5 * V_A(k, l, c, j) * t1(i, c) * t1(k, b) * t1(l, a);  // P(ab) applied
                    value += 0.5 * V_A(k, l, c, i) * t1(j, c) * t1(k, b) * t1(l, a);  // P(ij) P(ab) applied
                }
            }
        }

        for (const auto& k : orbital_space.occupiedIndices()) {
            for (const auto& c : orbital_space.virtualIndices()) {
                for (const auto& d : orbital_space.virtualIndices()) {
                    value += V_A(k, a, c, d) * t1(k, c) * t2(i, j, d, b);
                    value -= V_A(k, b, c, d) * t1(k, c) * t2(i, j, d, a);  // P(ab) applied

                    value += V_A(a, k, d, c) * t1(i, d) * t2(j, k, b, c);
                    value -= V_A(a, k, d, c) * t1(j, d) * t2(i, k, b, c);  // P(ij) applied
                    value -= V_A(b, k, d, c) * t1(i, d) * t2(j, k, a, c);  // P(ab) applied
                    value += V_A(b, k, d, c) * t1(j, d) * t2(i, k, a, c);  // P(ij) P(ab) applied

                    value -= 0.5 * V_A(k, b, c, d) * t1(k, a) * t2(i, j, c, d);
                    value += 0.5 * V_A(k, a, c, d) * t1(k, b) * t2(i, j, c, d);  // P(ab) applied

                    value -= 0.5 * V_A(k, b, c, d) * t1(i, c) * t1(k, a) * t1(j, d);
                    value += 0.5 * V_A(k, b, c, d) * t1(j, c) * t1(k, a) * t1(i, d);  // P(ij) applied
                    value += 0.5 * V_A(k, a, c, d) * t1(i, c) * t1(k, b) * t1(j, d);  // P(ab) applied
                    value -= 0.5 * V_A(k, a, c, d) * t1(j, c) * t1(k, b) * t1(i, d);  // P(ij) P(ab) applied
                }
            }
        }

        for (const auto& k : orbital_space.occupiedIndices()) {
            for (const auto& l : orbital_space.occupiedIndices()) {
                for (const auto& c : orbital_space.virtualIndices()) {
                    for (const auto& d : orbital_space.virtualIndices()) {
                        value += 0.5 * V_A(k, l, c, d) * t2(i, k, a, c) * t2(l, j, d, b);
                        value -= 0.5 * V_A(k, l, c, d) * t2(j, k, a, c) * t2(l, i, d, b);  // P(ij) applied
                        value -= 0.5 * V_A(k, l, c, d) * t2(i, k, b, c) * t2(l, j, d, a);  // P(ab) applied
                        value += 0.5 * V_A(k, l, c, d) * t2(j, k, b, c) * t2(l, i, d, a);  // P(ij) P(ab) applied

                        value += 0.25 * V_A(k, l, c, d) * t2(i, j, c, d) * t2(k, l, a, b);

                        value -= 0.5 * V_A(k, l, c, d) * t2(i, j, a, c) * t2(k, l, b, d);
                        value += 0.5 * V_A(k, l, c, d) * t2(i, j, b, c) * t2(k, l, a, d);  // P(ab) applied

                        value -= 0.5 * V_A(k, l, c, d) * t2(i, k, a, b) * t2(j, l, c, d);
                        value += 0.5 * V_A(k, l, c, d) * t2(j, k, a, b) * t2(i, l, c, d);  // P(ij) applied

                        value -= V_A(k, l, c, d) * t1(k, c) * t1(i, d) * t2(l, j, a, b);
                        value += V_A(k, l, c, d) * t1(k, c) * t1(j, d) * t2(l, i, a, b);  // P(ij) applied

                        value -= V_A(k, l, c, d) * t1(k, c) * t1(l, a) * t2(i, j, d, b);
                        value += V_A(k, l, c, d) * t1(k, c) * t1(l, b) * t2(i, j, d, a);  // P(ab) applied

                        value += 0.25 * V_A(k, l, c, d) * t1(i, c) * t1(j, d) * t2(k, l, a, b);
                        value -= 0.25 * V_A(k, l, c, d) * t1(j, c) * t1(i, d) * t2(k, l, a, b);  // P(ij) applied

                        value += 0.25 * V_A(k, l, c, d) * t1(k, a) * t1(l, b) * t2(i, j, c, d);
                        value -= 0.25 * V_A(k, l, c, d) * t1(k, b) * t1(l, a) * t2(i, j, c, d);  // P(ab) applied

                        value += V_A(k, l, c, d) * t1(i, c) * t1(l, b) * t2(k, j, a, d);
                        value -= V_A(k, l, c, d) * t1(j, c) * t1(l, b) * t2(k, i, a, d);  // P(ij) applied
                        value -= V_A(k, l, c, d) * t1(i, c) * t1(l, a) * t2(k, j, b, d);  // P(ab) applied
                        value += V_A(k, l, c, d) * t1(j, c) * t1(l, a) * t2(k, i, b, d);  // P(ij) P(ab) applied

                        value += 0.25 * V_A(k, l, c, d) * t1(i, c) * t1(k, a) * t1(j, d) * t1(l, b);
                        value -= 0.25 * V_A(k, l, c, d) * t1(j, c) * t1(k, a) * t1(i, d) * t1(l, b);  // P(ij) applied
                        value -= 0.25 * V_A(k, l, c, d) * t1(i, c) * t1(k, b) * t1(j, d) * t1(l, a);  // P(ab) applied
                        value += 0.25 * V_A(k, l, c, d) * t1(j, c) * t1(k, b) * t1(i, d) * t1(l, a);  // P(ij) P(ab) applied
                    }
                }
            }
        }


        return value;
    }
};


}  // namespace QCModel
}  // namespace GQCP
