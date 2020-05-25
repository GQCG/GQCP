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
#include "QCModel/CC/T2Amplitudes.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The CCD (coupled-cluster doubles) wave function model.
 * 
 *  @tparam _Scalar             the scalar type of the amplitudes
 */
template <typename _Scalar>
class CCD {
public:
    using Scalar = _Scalar;


private:
    T2Amplitudes<Scalar> t2;  // the T2-amplitudes


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Construct a CCD wave function from its converged T2-amplitudes.
     * 
     *  @param t2                   the T2-amplitudes
     */
    CCD(const T2Amplitudes<Scalar>& t2) :
        t2 {t2} {}


    /*
     *  STATIC PUBLIC METHODS
     */

    /**
     *  Calculate the CCD correlation energy.
     * 
     *  @param f                        the (inactive) Fock matrix
     *  @param V_A                      the antisymmetrized two-electron integrals (in physicist's notation)
     *  @param t2                       the T2-amplitudes
     * 
     *  @return the CCD correlation energy
     */
    static Scalar calculateCorrelationEnergy(const QCMatrix<Scalar>& f, const QCRankFourTensor<Scalar>& V_A, const T2Amplitudes<Scalar>& t2) {

        const auto orbital_space = t2.orbitalSpace(); 

        // A KISS implementation of the CCD energy correction.
        // The implementation is in line with Crawford2000 "Chapter 2: An Introduction to Coupled Cluster Theory for Computational Chemists", eq. [134].
        double E = 0.0;

        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                        E += 0.25 * V_A(i, j, a, b) * t2(i, j, a, b);
                    }
                }
            }
        }

        return E;
    }


    /**
     *  Calculate the value for one of the CCD T2-amplitude equations, evaluated at the given T2-amplitudes (and itermediates).
     *      f_{ij}^{ab} = <Phi_{ij}^{ab}| H |Phi_0>             with H the similarity-transformed normal-ordered Hamiltonian
     * 
     *  @param i                            the (occupied) subscript for the amplitude equation
     *  @param j                            the other (occupied) subscript for the amplitude equation
     *  @param a                            the (virtual) superscript for the amplitude equation
     *  @param b                            the other (virtual) superscript for the amplitude equation
     *  @param f                            the (inactive) Fock matrix
     *  @param V_A                          the antisymmetrized two-electron integrals (in physicist's notation)
     *  @param t2                           the T2-amplitudes
     *  @param F1                           the F1-intermediate (equation (3) in Stanton1991)
     *  @param F2                           the F2-intermediate (equation (4) in Stanton1991)
     *  @param F3                           the F3-intermediate (equation (5) in Stantion1991)
     *  @param W1                           the W1-intermediate (equation (6) in Stanton1991)
     *  @param W2                           the W2-intermediate (equation (7) in Stanton1991)
     *  @param W3                           the W3-intermediate (equation (8) in Stantion1991)
     * 
     *  @return the value for one of the CCD T2-amplitude equations
     */
    static Scalar calculateT2AmplitudeEquation(const size_t i, const size_t j, const size_t a, const size_t b, const QCMatrix<Scalar>& f, const QCRankFourTensor<Scalar>& V_A, const T2Amplitudes<Scalar>& t2, const ImplicitMatrixSlice<Scalar>& F1, const ImplicitMatrixSlice<Scalar>& F2, const ImplicitMatrixSlice<Scalar>& F3, const ImplicitRankFourTensorSlice<Scalar>& W1, const ImplicitRankFourTensorSlice<Scalar>& W2, const ImplicitRankFourTensorSlice<Scalar>& W3) {

        const auto orbital_space = t2.orbitalSpace(); 

        // We will use equation (2) in Stanton1991 by putting the left-hand term (with the energy denominator) to the right.
        Scalar result {0.0};  // zero-initialize the scalar value for the result

        // Calculate the contribution from the left-hand side.
        const auto D_ijab = f(i, i) + f(j, j) - f(a, a) - f(b, b);
        result -= t2(i, j, a, b) * D_ijab;

        // Calculate the contribution from the first term.
        result += V_A(i, j, a, b);

        // Calculate the contribution from the second term.
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
            result += t2(i, j, a, e) * F1(b, e);
            result -= t2(i, j, b, e) * F1(a, e);  // P(ab) applied
        }

        // Calculate the contribution from the third term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            result -= t2(i, m, a, b) * F2(m, j);
            result += t2(j, m, a, b) * F2(m, i);  // P(ij) applied
        }

        // Calculate the contribution from the fourth term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                result += 0.5 * t2(m, n, a, b) * W1(m, n, i, j);
            }
        }

        // Calculate the contribution from the fifth term.
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                result += 0.5 * t2(i, j, e, f) * W2(a, b, e, f);
            }
        }

        // Calculate the contribution from the sixth term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                result += t2(i, m, a, e) * W3(m, b, e, j);
                result -= t2(j, m, a, e) * W3(m, b, e, i);  // P(ij) applied
                result -= t2(i, m, b, e) * W3(m, a, e, j);  // P(ab) applied
                result += t2(j, m, b, e) * W3(m, a, e, i);  // P(ij) P(ab) applied
            }
        }

        return result;
    }


    /**
     *  @param f                    the (inactive) Fock matrix
     *  @param V_A                  the antisymmetrized two-electron integrals (in physicist's notation)
     *  @param t2                   the T2-amplitudes
     * 
     *  @return the F1-intermediate
     * 
     *  @note This is one of the intermediate quantities in the factorization of CCD. In particular, F1 represents equation (3) in Stanton1991.
     */
    static ImplicitMatrixSlice<Scalar> calculateF1(const QCMatrix<Scalar>& f, const QCRankFourTensor<Scalar>& V_A, const T2Amplitudes<Scalar>& t2) {

        const auto& orbital_space = t2.orbitalSpace();

        // Implement the formula for the F1-intermediate: equation (3) in Stanton1993.
        auto F1 = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_virtual, OccupationType::k_virtual);  // zero-initialize a virtual-virtual object
        for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                Scalar value {0.0};  // zero-initialize the scalar value to be added

                // Calculate the contribution from the first term.
                if (a != e) {  // (1 - delta_ae)
                    value += f(a, e);
                }

                // Calculate the contribution from the fourth term.
                for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                        for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                            value -= 0.5 * t2(m, n, a, f) * V_A(m, n, e, f);
                        }
                    }
                }

                F1(a, e) = value;
            }
        }

        return F1;
    }


    /**
     *  @param f                    the (inactive) Fock matrix
     *  @param V_A                  the antisymmetrized two-electron integrals (in physicist's notation)
     *  @param t2                   the T2-amplitudes
     * 
     *  @return the F2-intermediate
     * 
     *  @note This is one of the intermediate quantities in the factorization of CCD. In particular, F2 represents equation (4) in Stanton1991.
     */
    static ImplicitMatrixSlice<Scalar> calculateF2(const QCMatrix<Scalar>& f, const QCRankFourTensor<Scalar>& V_A, const T2Amplitudes<Scalar>& t2) {

        const auto& orbital_space = t2.orbitalSpace();

        // Implement the formula for F2 in equation (4) in Stanton1991.
        auto F2 = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_occupied);  // zero-initialize an occupied-occupied object
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
                Scalar value {0.0};  // zero-initialize the scalar value to be added

                // Calculate the contribution from the first term.
                if (m != i) {  // (1 - delta_mi)
                    value += f(m, i);
                }

                // Calculate the contribution from the fourth term.
                for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                        for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                            value += 0.5 * t2(i, n, e, f) * V_A(m, n, e, f);
                        }
                    }
                }

                F2(m, i) = value;
            }
        }

        return F2;
    }


    /**
     *  @param f                    the (inactive) Fock matrix
     *  @param V_A                  the antisymmetrized two-electron integrals (in physicist's notation)
     *  @param t2                   the T2-amplitudes
     * 
     *  @return the F3-intermediate
     * 
     *  @note This is one of the intermediate quantities in the factorization of CCD. In particular, F3 represents equation (5) in Stanton1991.
     */
    static ImplicitMatrixSlice<Scalar> calculateF3(const QCMatrix<Scalar>& f, const QCRankFourTensor<Scalar>& V_A, const T2Amplitudes<Scalar>& t2) {

        const auto& orbital_space = t2.orbitalSpace();

        // Implement the formula for F3 in equation (5) in Stanton1991.
        auto F3 = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual);  // zero-initialize an occupied-virtual object
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                // Calculate the contribution from the first term.
                F3(m, e) = f(m, e);
            }
        }

        return F3;
    }


    /**
     *  @param V_A                  the antisymmetrized two-electron integrals (in physicist's notation)
     *  @param t2                   the T2-amplitudes
     * 
     *  @return the W1-intermediate
     * 
     *  @note This is one of the intermediate quantities in the factorization of CCD. In particular, W1 represents equation (6) in Stanton1991.
     */
    static ImplicitRankFourTensorSlice<Scalar> calculateW1(const QCRankFourTensor<Scalar>& V_A, const T2Amplitudes<Scalar>& t2) {

        const auto& orbital_space = t2.orbitalSpace();

        // Implement the formula for W1 (equation 6).
        auto W1 = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_occupied, OccupationType::k_occupied, OccupationType::k_occupied);  // zero-initialize an occupied-occupied-occupied-occupied object
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto j : orbital_space.indices(OccupationType::k_occupied)) {
                        Scalar value {0.0};  // zero-initialize the scalar value to be added

                        // Calculate the contribution from the first term.
                        value += V_A(m, n, i, j);

                        // Calculate the contribution from the third term.
                        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                            for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                                value += 0.25 * t2(i, j, e, f) * V_A(m, n, e, f);
                            }
                        }

                        W1(m, n, i, j) = value;
                    }
                }
            }
        }

        return W1;
    }


    /**
     *  @param V_A                  the antisymmetrized two-electron integrals (in physicist's notation)
     *  @param t2                   the T2-amplitudes
     * 
     *  @return the W2-intermediate
     * 
     *  @note This is one of the intermediate quantities in the factorization of CCD. In particular, W2 represents equation (7) in Stanton1991.
     */
    static ImplicitRankFourTensorSlice<Scalar> calculateW2(const QCRankFourTensor<Scalar>& V_A, const T2Amplitudes<Scalar>& t2) {

        const auto& orbital_space = t2.orbitalSpace();

        // Implement the formula for W2 (equation 7).
        auto W2 = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_virtual, OccupationType::k_virtual, OccupationType::k_virtual, OccupationType::k_virtual);  // zero-initialize a virtual-virtual-virtual-virtual object
        for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto f : orbital_space.indices(OccupationType::k_virtual)) {
                        Scalar value {0.0};  // zero-initialize the scalar value to be added

                        // Calculate the contribution from the first term.
                        value += V_A(a, b, e, f);

                        // Calculate the contribution from the third term.
                        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
                            for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                                value += 0.25 * t2(m, n, a, b) * V_A(m, n, e, f);
                            }
                        }

                        W2(a, b, e, f) = value;
                    }
                }
            }
        }

        return W2;
    }


    /**
     *  @param V_A                  the antisymmetrized two-electron integrals (in physicist's notation)
     *  @param t2                   the T2-amplitudes
     * 
     *  @return the W3-intermediate
     * 
     *  @note This is one of the intermediate quantities in the factorization of CCD. In particular, W3 represents equation (8) in Stanton1991.
     */
    static ImplicitRankFourTensorSlice<Scalar> calculateW3(const QCRankFourTensor<Scalar>& V_A, const T2Amplitudes<Scalar>& t2) {

        const auto& orbital_space = t2.orbitalSpace();

        // Implement the formula for W3 (equation 8).
        auto W3 = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_virtual, OccupationType::k_occupied);  // zero-initialize an occupied-virtual-virtual-occupied object
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                        Scalar value {0.0};  // zero-initialize the scalar value to be added

                        // Calculate the contribution from the first term.
                        value += V_A(m, b, e, j);

                        // Calculate the contribution from the fourth term.
                        for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                            for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                                value -= 0.5 * t2(j, n, f, b);
                            }
                        }

                        W3(m, b, e, j) = value;
                    }
                }
            }
        }

        return W3;
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return these CCD model parameters' T2-amplitudes
     */
    const T2Amplitudes<Scalar>& t2Amplitudes() const { return this->t2; }
};


}  // namespace QCModel
}  // namespace GQCP
