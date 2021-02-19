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
#include "QCModel/CC/L1Amplitudes.hpp"
#include "QCModel/CC/L2Amplitudes.hpp"
#include "QCModel/CC/T1Amplitudes.hpp"
#include "QCModel/CC/T2Amplitudes.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The LambdaCCSD (Lambda-equation coupled-cluster singles and doubles) wave function model.
 * 
 *  @tparam _Scalar             The scalar type of the amplitudes.
 */
template <typename _Scalar>
class LambdaCCSD {
public:
    // The scalar type of the amplitudes.
    using Scalar = _Scalar;


private:
    // The T1-amplitudes.
    T1Amplitudes<Scalar> t1;

    // The L1-amplitudes
    L1Amplitudes<Scalar> l1;

    // The T2-amplitudes.
    T2Amplitudes<Scalar> t2;

    // The L2-amplitudes.
    L2Amplitudes<Scalar> l2;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a LambdaCCSD wave function from its converged T1- and T2-amplitudes.
     * 
     *  @param t1                   The T1-amplitudes.
     *  @param t2                   The T2-amplitudes.
     *  @param l1                   The L1-amplitudes.
     *  @param l2                   The L2-amplitudes.
     */
    LambdaCCSD(const T1Amplitudes<Scalar>& t1, const T2Amplitudes<Scalar>& t2, const L1Amplitudes<Scalar>& l1, const L2Amplitudes<Scalar>& l2) :
        t1 {t1},
        l1 {l1},
        t2 {t2},
        l2 {l2} {}


    /*
     *  MARK: Access
     */

    /**
     *  @return These LambdaCCSD model parameters' T1-amplitudes.
     */
    const T1Amplitudes<Scalar>& t1Amplitudes() const { return this->t1; }

    /**
     *  @return These LambdaCCSD model parameters' T2-amplitudes.
     */
    const T2Amplitudes<Scalar>& t2Amplitudes() const { return this->t2; }

    /**
     *  @return These LambdaCCSD model parameters' L1-amplitudes.
     */
    const T1Amplitudes<Scalar>& l1Amplitudes() const { return this->l1; }

    /**
     *  @return These LambdaCCSD model parameters' L2-amplitudes.
     */
    const T2Amplitudes<Scalar>& l2Amplitudes() const { return this->l2; }


    /*
     *  MARK: Intermediates for Lambda-CCSD
     */

    /**
     *  @param F1           The F1-intermediate (equation (4) in Gauss1991).
     *  @param F3           The F3-intermediate (equation (6) in Gauss1991).
     *  @param t1           The T1-amplitudes.
     * 
     *  @return the F1_tilde-intermediate
     * 
     *  @note This is one of the intermediate quantities in the factorization of Lambda-CCSD. In particular, F1_tilde represents equation (17) in Gauss1991.
     */
    static ImplicitMatrixSlice<Scalar> calculateF1tilde(const ImplicitMatrixSlice<Scalar>& F1, const ImplicitMatrixSlice<Scalar>& F3, const T1Amplitudes<Scalar>& t1) {

        const auto& orbital_space = t1.orbitalSpace();

        auto F1_tilde = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_virtual, OccupationType::k_virtual);  // Zero-initialize a virtual-virtual object.
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                // Calculate the contribution from the first term.
                Scalar value = F1(e, a);

                // Calculate the contribution from the second term.
                for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
                    value -= 0.5 * t1(m, e) * F3(m, a);
                }

                F1_tilde(e, a) = value;
            }
        }

        return F1_tilde;
    }


    /**
     *  @param F2           The F2-intermediate (equation (5) in Gauss1991).
     *  @param F3           The F3-intermediate (equation (6) in Gauss1991).
     *  @param t1           The T1-amplitudes.
     * 
     *  @return The F2_tilde-intermediate.
     * 
     *  @note This is one of the intermediate quantities in the factorization of Lambda-CCSD. In particular, F2_tilde represents equation (18) in Gauss1991.
     */
    static ImplicitMatrixSlice<Scalar> calculateF2tilde(const ImplicitMatrixSlice<Scalar>& F2, const ImplicitMatrixSlice<Scalar>& F3, const T1Amplitudes<Scalar>& t1) {

        const auto& orbital_space = t1.orbitalSpace();

        auto F2_tilde = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_occupied);  // Zero-initialize a virtual-virtual object.
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
                // Calculate the contribution from the first term.
                Scalar value = F2(i, m);

                // Calculate the contribution from the second term.
                for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                    value -= 0.5 * t1(m, e) * F3(i, e);
                }

                F2_tilde(i, m) = value;
            }
        }

        return F2_tilde;
    }


    /**
     *  @param V_A                  The antisymmetrized two-electron integrals (in physicist's notation).
     *  @param tau2                 The tau2-intermediate (equation (11) in Gauss1991).
     *  @param W2                   The W2-intermediate (equation (8) in Gauss1991).
     *  @param orbital_space        The orbital space that encapsulates the occupied-virtual separation.
     * 
     *  @return The W1_tilde-intermediate.
     * 
     *  @note This is one of the intermediate quantities in the factorization of Lambda-CCSD. In particular, W1_tilde represents equation (19) in Gauss1991.
     */
    static ImplicitRankFourTensorSlice<Scalar> calculateW1Tilde(const SquareRankFourTensor<Scalar>& V_A, const ImplicitRankFourTensorSlice<Scalar>& tau2, const ImplicitRankFourTensorSlice<Scalar>& W2, const OrbitalSpace& orbital_space) {

        auto W1_tilde = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_virtual, OccupationType::k_virtual, OccupationType::k_virtual, OccupationType::k_virtual);  // Zero-initialize a virtual-virtual-virtual-virtual object.
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {

                        // Calculate the contribution from the first term.
                        Scalar value = W2(e, f, a, b);

                        // Calculate the contribution from the second term.
                        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
                            for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                                value += 0.25 * tau2(m, n, e, f) * V_A(m, n, a, b);
                            }
                        }

                        W1_tilde(a, b, e, f) = value;
                    }
                }
            }
        }

        return W1_tilde;
    }


    /**
     *  @param V_A                  The antisymmetrized two-electron integrals (in physicist's notation).
     *  @param tau2                 The tau2-intermediate (equation (11) in Gauss1991).
     *  @param W2                   The W2-intermediate (equation (8) in Gauss1991).
     *  @param orbital_space        The orbital space that encapsulates the occupied-virtual separation.
     * 
     *  @return The W2_tilde-intermediate.
     * 
     *  @note This is one of the intermediate quantities in the factorization of Lambda-CCSD. In particular, W2_tilde represents equation (20) in Gauss1991.
     */
    static ImplicitRankFourTensorSlice<Scalar> calculateW2Tilde(const SquareRankFourTensor<Scalar>& V_A, const ImplicitRankFourTensorSlice<Scalar>& tau2, const ImplicitRankFourTensorSlice<Scalar>& W2, const OrbitalSpace& orbital_space) {

        auto W2_tilde = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_occupied, OccupationType::k_occupied, OccupationType::k_occupied);  // Zero-initialize an occupied-occupied-occupied-occupied object.
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {

                        // Calculate the contribution from the first term.
                        Scalar value = W1(i, j, m, n);

                        // Calculate the contribution from the second term.
                        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                            for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                                value += 0.25 * tau2(m, n, e, f) * V_A(e, f, i, j);
                            }
                        }

                        W2_tilde(i, j, m, n) = value;
                    }
                }
            }
        }

        return W2_tilde;
    }


    /**
     *  @param V_A                  The antisymmetrized two-electron integrals (in physicist's notation).
     *  @param t1                   The T1-amplitudes.
     *  @param t2                   The T2-amplitudes.
     * 
     *  @return The W3_tilde-intermediate.
     * 
     *  @note This is one of the intermediate quantities in the factorization of Lambda-CCSD. In particular, W3_tilde represents equation (21) in Gauss1991.
     */
    static ImplicitRankFourTensorSlice<Scalar> calculateW3Tilde(const SquareRankFourTensor<Scalar>& V_A, const T1Amplitudes<Scalar>& t1, const T2Amplitudes<Scalar>& t2) {

        const auto& orbital_space = t1.orbitalSpace();  // Assume the orbital spaces for t1 and t2 are equal.

        auto W3_tilde = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_occupied, OccupationType::k_virtual);  // Zero-initialize a virtual-occupied-occupied-virtual object
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {

                        // Zero-initialize the scalar value to be added.
                        Scalar value {0.0};

                        // Calculate the contribution from the first term.
                        value += V_A(e, j, m, b);

                        // Calculate the contribution from the second term.
                        for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                            value += t1(m, f) * V_A(e, j, f, b);
                        }

                        // Calculate the contribution from the third term.
                        for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                            value -= t1(n, e) * V_A(n, j, m, b);
                        }

                        // Calculate the contribution from the fourth term.
                        for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                            for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                                value += (t2(m, n, e, f) - t1(m, f) * t1(n, e)) * V_A(n, j, f, b);
                            }
                        }

                        W3_tilde(e, j, m, b) = value;
                    }
                }
            }
        }

        return W3_tilde;
    }


    /**
     *  @param V_A                  The antisymmetrized two-electron integrals (in physicist's notation).
     *  @param t1                   The T1-amplitudes.
     *  @param t2                   The T2-amplitudes.
     *  @param tau2
     *  @param F3
     *  @param W3_tilde
     *  @param W_tilde_tilde
     * 
     *  @return The W4_tilde-intermediate.
     * 
     *  @note This is one of the intermediate quantities in the factorization of Lambda-CCSD. In particular, W4_tilde represents equation (22) in Gauss1991.
     */
    static ImplicitRankFourTensorSlice<Scalar> calculateW4Tilde(const SquareRankFourTensor<Scalar>& V_A, const T1Amplitudes<Scalar>& t1, const T2Amplitudes<Scalar>& t2, const ImplicitRankFourTensorSlice<Scalar>& tau2, const ImplicitMatrixSlice<Scalar>& F3, const ImplicitRankFourTensorSlice<Scalar>& W3_tilde, const ImplicitRankFourTensorSlice<Scalar>& W_tilde_tilde) {

        const auto& orbital_space = t1.orbitalSpace();  // Assume the orbital spaces for t1 and t2 are equal.

        auto W4_tilde = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_occupied);  // Zero-initialize a virtual-occupied-occupied-virtual object
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {

                        // Zero-initialize the scalar value to be added.
                        Scalar value {0.0};

                        // Calculate the contribution from the first term.
                        value += V_A(i, e, m, n);

                        // Calculate the contribution from the second term.
                        for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                            value -= F3(i, f) * t2(m, n, e, f);
                        }

                        // Calculate the contribution from the third term.
                        for (const auto& o : orbital_space.indices(OccupationType::k_occupied)) {
                            value -= t1(o, e) * W3_tilde(i, o, m, n);
                        }

                        // Calculate the contribution from the fourth term.
                        for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                            for (const auto& g : orbital_space.indices(OccupationType::k_virtual)) {
                                value += 0.5 * V_A(i, e, f, g) * tau2(m, n, f, g);
                            }
                        }

                        // Calculate the contribution from the fifth term.
                        for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                            value -= t1(m, f) * W_tilde_tilde(i, e, f, n);
                            value += t1(n, f) * W_tilde_tilde(i, e, f, m);  // P(mn) applied.
                        }

                        // Calculate the contribution from the sixth term.
                        for (const auto& o : orbital_space.indices(OccupationType::k_occupied)) {
                            for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                                value += V_A(i, o, m, f) * t2(n, o, e, f);
                                value -= V_A(i, o, n, f) * t2(m, o, e, f);  // P(mn) applied.
                            }
                        }

                        W4_tilde(e, j, m, b) = value;
                    }
                }
            }
        }

        return W4_tilde;
    }


    /**
     *  @param V_A                  The antisymmetrized two-electron integrals (in physicist's notation).
     *  @param t1                   The T1-amplitudes.
     *  @param t2                   The T2-amplitudes.
     *  @param tau2
     *  @param F3
     *  @param W1_tilde
     *  @param W_tilde_tilde
     * 
     *  @return The W5_tilde-intermediate.
     * 
     *  @note This is one of the intermediate quantities in the factorization of Lambda-CCSD. In particular, W5_tilde represents equation (23) in Gauss1991.
     */
    static ImplicitRankFourTensorSlice<Scalar> calculateW5Tilde(const SquareRankFourTensor<Scalar>& V_A, const T1Amplitudes<Scalar>& t1, const T2Amplitudes<Scalar>& t2, const ImplicitRankFourTensorSlice<Scalar>& tau2, const ImplicitMatrixSlice<Scalar>& F3, const ImplicitRankFourTensorSlice<Scalar>& W1_tilde, const ImplicitRankFourTensorSlice<Scalar>& W_tilde_tilde) {

        const auto& orbital_space = t1.orbitalSpace();  // Assume the orbital spaces for t1 and t2 are equal.

        auto W5_tilde = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_virtual, OccupationType::k_virtual, OccupationType::k_virtual, OccupationType::k_occupied);  // Zero-initialize a virtual-virtual-virtual-occupied object
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {

                        // Zero-initialize the scalar value to be added.
                        Scalar value {0.0};

                        // Calculate the contribution from the first term.
                        value += V_A(e, f, a, m);

                        // Calculate the contribution from the second term.
                        for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                            value += t2(m, n, e, f) * F3(n, a);
                        }

                        // Calculate the contribution from the third term.
                        for (const auto& g : orbital_space.indices(OccupationType::k_virtual)) {
                            value += t1(m, g) * W1_tilde(e, f, a, g);
                        }

                        // Calculate the contribution from the fourth term.
                        for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                            for (const auto& o : orbital_space.indices(OccupationType::k_occupied)) {
                                value += 0.5 * V_A(a, m, n, o) * tau2(n, o, e, f);
                            }
                        }

                        // Calculate the contribution from the fifth term.
                        for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                            value -= t1(n, e) * W_tilde_tilde(n, f, a, m);
                            value += t1(n, f) * W_tilde_tilde(n, e, a, m);  // P(ef) applied.
                        }

                        // Calculate the contribution from the sixth term.
                        for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                            for (const auto& g : orbital_space.indices(OccupationType::k_virtual)) {
                                value += V_A(e, n, a, g) * t2(m, n, f, g);
                                value -= V_A(f, n, a, g) * t2(m, n, e, g);  // P(ef) applied
                            }
                        }

                        W5_tilde(e, j, m, b) = value;
                    }
                }
            }
        }

        return W5_tilde;
    }


    /**
     *  @param V_A                  The antisymmetrized two-electron integrals (in physicist's notation).
     *  @param t2                   The T2-amplitudes.
     * 
     *  @return The W_tilde_tilde-intermediate.
     * 
     *  @note This is one of the intermediate quantities in the factorization of Lambda-CCSD. In particular, W_tilde_tilde represents equation (24) in Gauss1991.
     */
    static ImplicitRankFourTensorSlice<Scalar> calculateWTildeTilde(const SquareRankFourTensor<Scalar>& V_Aconst T2Amplitudes<Scalar>& t2) {

        const auto& orbital_space = t2.orbitalSpace();

        auto W_tilde_tilde = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_virtual, OccupationType::k_occupied);  // Zero-initialize a occupied-virtual-virtual-occupied object.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {

                        // Calculate the contribution from the first term.
                        Scalar value = V_A(m, b, e, f);

                        // Calculate the contribution from the second term.
                        for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                            for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                                value -= t2(n, j, b, f) * V_A(m, n, e, f);
                            }
                        }

                        W_tilde_tilde(m, b, e, j) = value;
                    }
                }
            }
        }

        return W_tilde_tilde;
    }


    /**
     *  @param t2           The T2-amplitudes.
     *  @param l2           The L2-amplitudes.
     * 
     *  @return The G1-intermediate.
     * 
     *  @note This is one of the intermediate quantities in the factorization of Lambda-CCSD. In particular, G1 represents equation (25) in Gauss1991.
     */
    static ImplicitMatrixSlice<Scalar> calculateG1(const T2Amplitudes<Scalar>& t2, const L2Amplitudes<Scalar>& l2) {

        const auto& orbital_space = t2.orbitalSpace();  // Assume the orbital spaces for t2 and l2 are equal.

        auto G1 = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_virtual, OccupationType::k_virtual);  // Zero-initialize a virtual-virtual object.
        for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {

                Scalar value {0.0};

                for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                        for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                            value -= 0.5 * t2(m, n, e, f) * l2(a, f, m, n);
                        }
                    }
                }

                G1(a, e) = value;
            }
        }

        return G1;
    }


    /**
     *  @param t2           The T2-amplitudes.
     *  @param l2           The L2-amplitudes.
     * 
     *  @return The G2-intermediate.
     * 
     *  @note This is one of the intermediate quantities in the factorization of Lambda-CCSD. In particular, G2 represents equation (26) in Gauss1991.
     */
    static ImplicitMatrixSlice<Scalar> calculateG1(const T2Amplitudes<Scalar>& t2, const L2Amplitudes<Scalar>& l2) {

        const auto& orbital_space = t2.orbitalSpace();  // Assume the orbital spaces for t2 and l2 are equal.

        auto G2 = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_occupied);  // Zero-initialize an occupied-occupied object.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {

                Scalar value {0.0};

                for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                        for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                            value += t2(m, n, e, f) * l2(e, f, i, n);
                        }
                    }
                }

                G2(m, i) = value;
            }
        }

        return G2;
    }


    /*
     *  MARK: Lambda equations
     */

    /**
     *  Calculate the value for one of the Lambda-CCSD L1-amplitude equations, evaluated at the intermediates.
     * 
     *  @param a                    The virtual index for the amplitude equation.
     *  @param i                    The occupied index for the amplitude equation.
     *  @param f                    The (inactive) Fock matrix.
     *  @param V_A                  The antisymmetrized two-electron integrals (in physicist's notation).
     * 
     * 
     * 
     * 
     * 
     * 
     * 
     * 
     *  @return The value for one of the CCSD L1-amplitude equations.
     */
    static Scalar calculateT1AmplitudeEquation(const size_t a, const size_t i, const SquareMatrix<Scalar>& f, const SquareRankFourTensor<Scalar>& V_A, const T1Amplitudes<Scalar>& t1, const L1Amplitudes<Scalar>& l1, const L2Amplitudes<Scalar>& l2, const ImplicitMatrixSlice<Scalar>& G1, const ImplicitMatrixSlice<Scalar>& G2, const ImplicitMatrixSlice<Scalar>& F1_tilde, const ImplicitMatrixSlice<Scalar>& F2_tilde, const ImplicitRankFourTensorSlice<Scalar>& W3_tilde, const ImplicitRankFourTensorSlice<Scalar>& W4_tilde, const ImplicitRankFourTensorSlice<Scalar>& W5_tilde) {

        const auto& orbital_space = t1.orbitalSpace();  // Assume all implements have the same orbital space.

        // We will use equation (15) in Gauss1991 by putting the left-hand term (with the energy denominator) to the right.
        Scalar result {0.0};  // Zero-initialize the scalar value for the result.

        // Calculate the contribution from the left-hand side.
        const auto D_ia = f(i, i) - f(a, a);
        result -= D_ia * l1(a, i);

        // Calculate the contribution from the first term.
        value += F3(i, a);

        // Calculate the contribution from the second term.
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
            value += l1(e, i) * F1_tilde(e, a);
        }

        // Calculate the contribution from the third term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            value -= l1(a, m) * F2_tilde(i, m);
        }

        // Calculate the contribution from the fourth term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                    value += 0.5 * l2(e, f, i, m) * W5_tilde(e, f, a, m);
                }
            }
        }

        // Calculate the contribution from the fifth term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                value += l1(e, m) * W3_tilde(e, i, m, a);
            }
        }

        // Calculate the contribution from the sixth term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                    value -= 0.5 * l2(a, e, m, n) * W4_tilde(i, e, m, n);  // We believe the order of the indices that is provided in the paper, i.e. W_mnie, is wrong and should be permuted.
                }
            }
        }

        // Calculate the contribution from the seventh term.
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                value -= G1(e, f) * V_A(e, i, f, a);
            }
        }

        // Calculate the contribution from the eight term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                value -= G2(m, n) * V_A(m, i, n, a);
            }
        }

        // Calculate the contribution from the ninth term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {

                for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                    value += G1(f, e) * t1(m, f) * V_A(i, m, a, e);
                }

                for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                    value -= G2(m, n) * t1(m, e) * V_A(i, m, a, e);
                }
            }
        }

        return result;
    }


    /**
     *  Calculate the value for one of the Lambda-CCSD L2-amplitude equations, evaluated at the intermediates.
     * 
     *  @param a                    The virtual index for the amplitude equation.
     *  @param i                    The occupied index for the amplitude equation.
     *  @param f                    The (inactive) Fock matrix.
     *  @param V_A                  The antisymmetrized two-electron integrals (in physicist's notation).
     * 
     * 
     * 
     * 
     * 
     * 
     * 
     * 
     *  @return The value for one of the CCSD L2-amplitude equations.
     */
    static Scalar
    calculateL2AmplitudeEquation(const size_t a, const size_t b, const size_t i, const size_t j, const SquareMatrix<Scalar>& f, const SquareRankFourTensor<Scalar>& V_A, const T1Amplitudes<Scalar>& t1, const L1Amplitudes<Scalar>& l1, const L2Amplitudes<Scalar>& l2, const ImplicitMatrixSlice<Scalar>& G1, const ImplicitMatrixSlice<Scalar>& G2, const ImplicitMatrixSlice<Scalar>& F3, const ImplicitMatrixSlice<Scalar>& F1_tilde, const ImplicitMatrixSlice<Scalar>& F2_tilde, const ImplicitRankFourTensorSlice<Scalar>& W1_tilde, const ImplicitRankFourTensorSlice<Scalar>& W2_tilde, const ImplicitRankFourTensorSlice<Scalar>& W3_tilde) {

        const auto& orbital_space = t1.orbitalSpace();  // Assume all implements have the same orbital space.

        // We will use equation (16) in Gauss1991 by putting the left-hand term (with the energy denominator) to the right.
        Scalar result {0.0};  // Zero-initialize the scalar value for the result.

        // Calculate the contribution from the left-hand side.
        const auto D_ijab = f(i, i) + f(j, j) - f(a, a) - f(b, b);
        result -= D_ijab * l2(a, b, i, j);

        // Calculate the contribution from the first term.
        value += V_A(i, j, a, b);

        // Calculate the contribution from the second term.
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
            // We believe the order of the indices that is provided in the paper, i.e. l2_ijae, is wrong and should be permuted.
            value += l2(a, e, i, j) * F1_tilde(e, b);
            value -= l2(b, e, i, j) * F1_tilde(e, a);  // P(ab) applied.
        }

        // Calculate the contribution from the third term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            value -= l2(a, b, i, m) * F2_tilde(j, m);
            value += l2(a, b, j, m) * F2_tilde(i, m);  // P(ij) applied.
        }

        // Calculate the contribution from the fourth term.
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& f : orbital_space.indices(OccupationType::k_virtual)) {
                value += 0.5 * l2(e, f, i, j) * W1_tilde(e, f, a, b);
            }
        }

        // Calculate the contribution from the fifth term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& n : orbital_space.indices(OccupationType::k_occupied)) {
                value += 0.5 * l2(a, b, m, n) * W2_tilde(i, j, m, n);
            }
        }

        // Calculate the contribution from the sixth term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                value += l2(a, e, i, m) * W3_tilde(e, j, m, b);
                value -= l2(a, e, j, m) * W3_tilde(e, i, m, b);  // P(ij) appplied.
                value -= l2(b, e, i, m) * W3_tilde(e, j, m, a);  // P(ab) applied.
                value += l2(b, e, j, m) * W3_tilde(e, i, m, a);  // P(ij) P(ab) applied.
            }
        }

        // Calculate the contribution from the seventh term.
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {

            // First contribution in parentheses.
            value += V_A(i, j, a, e) * G1(b, e);

            // Second contribution in parentheses.
            for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
                value -= V_A(i, j, a, e) * l1(b, m) * t1(m, e);
            }

            // First contribution in parentheses, with P(ab) applied.
            value -= V_A(i, j, b, e) * G1(a, e);

            // Second contribution in parentheses, with P(ab) applied.
            for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {
                value += V_A(i, j, b, e) * l1(a, m) * t1(m, e);
            }
        }

        // Calculate the contribution from the eight term.
        for (const auto& m : orbital_space.indices(OccupationType::k_occupied)) {

            // First contribution in parentheses.
            value -= V_A(i, m, a, b) * G2(m, j);

            // Second contribution in parentheses.
            for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                value -= V_A(i, m, a, b) * l1(e, j) * t1(m, e);
            }

            // First contribution in parentheses, with P(ij) applied.
            value += V_A(j, m, a, b) * G2(m, i);

            // Second contribution in parentheses, with P(ij) applied.
            for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
                value += V_A(j, m, a, b) * l1(e, i) * t1(m, e);
            }
        }

        // Calculate the contribution from the ninth term.
        value += l1(a, i) * F3(j, b);
        value -= l1(a, j) * F3(i, b);  // P(ij) applied
        value -= l1(b, i) * F3(j, a);  // P(ab) applied
        value += l1(b, j) * F3(i, a);  // P(ij) P(ab) applied.

        // Calculate the contribution from the eleventh term.
        for (const auto& e : orbital_space.indices(OccupationType::k_virtual)) {
            value += l1(e, i) * V_A(e, j, a, b);
            value -= l1(e, j) * V_A(e, i, a, b);  // P(ij) applied.
        }

        // Calculate the contribution from the twelfth term.
        for (const auto& e : orbital_space.indices(OccupationType::k_occupied)) {
            value -= l1(a, m) * V_A(i, j, m, b);
            value += l1(b, m) * V_A(i, j, m, a);  // P(ab) applied.
        }

        return value;
    }

    /*
     *  MARK: Density matrices
     */
};


}  // namespace QCModel
}  // namespace GQCP
