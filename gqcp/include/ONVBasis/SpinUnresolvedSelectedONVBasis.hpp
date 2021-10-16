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

#include "Mathematical/Representation/MatrixRepresentationEvaluationContainer.hpp"
#include "ONVBasis/SpinUnresolvedONV.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/GSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Utilities/complex.hpp"


namespace GQCP {


/**
 *  A spin-unresolved ONV basis with a flexible number of (spin-unresolved) ONVs.
 */
class SpinUnresolvedSelectedONVBasis {
public:
    // The ONV that is naturally related to a full spin-resolved ONV basis.
    using ONV = SpinUnresolvedONV;

private:
    // The number of spinors.
    size_t M;

    // The number of electrons.
    size_t N;

    // A collection of ONVs that span a 'selected' part of a Fock space.
    std::vector<SpinUnresolvedONV> onvs;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct an empty spin-unresolved selected ONV basis.
     *
     *  @param M            The number of spinors.
     *  @param N            The number of electrons.
     */
    SpinUnresolvedSelectedONVBasis(const size_t M, const size_t N);

    /**
     *  Generate a `SpinUnresolvedSelectedONVBasis` from a full spin-unresolved ONV basis.
     *
     *  @param onv_basis        The full spin-unresolved ONV basis.
     */
    SpinUnresolvedSelectedONVBasis(const SpinUnresolvedONVBasis& onv_basis);


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of spinors.
     */
    size_t numberOfOrbitals() const { return this->M; }

    /**
     *  @return The number of electrons.
     */
    size_t numberOfElectrons() const { return this->N; }

    /**
     *  @return The dimension of the Fock subspace that is spanned by this selected ONV basis.
     */
    size_t dimension() const { return this->onvs.size(); }


    /*
     *  MARK: Modifying
     */

    /**
     *  Expand this ONV basis with the given spin-unresolved ONV.
     *
     *  @param onv          The ONV that should be included in this ONV basis.
     */
    void expandWith(const SpinUnresolvedONV& onv);

    /**
     *  Expand this ONV basis with the given spin-unresolved ONVs.
     *
     *  @param onvs         The ONVs that should be included in this ONV basis.
     */
    void expandWith(const std::vector<SpinUnresolvedONV>& onvs);


    /*
     *  MARK: Accessing
     */

    /**
     *  Access the ONV that corresponds to the given index/address.
     *
     *  @param index            The address of the ONV.
     *
     *  @return The ONV that corresponds to the given index/address.
     */
    const SpinUnresolvedONV& onvWithIndex(const size_t index) const { return this->onvs[index]; }


    /*
     *  MARK: Dense generalized operator evaluations
     */

    /**
     *  Calculate the dense matrix representation of a generalized one-electron operator in this ONV basis.
     *
     *  @tparam Scalar          The scalar representation of a one-electron parameter: real or complex.
     *
     *  @param f                A generalized one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the one-electron operator.
     */
    template <typename Scalar>
    SquareMatrix<Scalar> evaluateOperatorDense(const ScalarGSQOneElectronOperator<Scalar>& f) const {

        if (f.numberOfOrbitals() != this->numberOfOrbitals()) {
            throw std::invalid_argument("SpinUnresolvedSelectedONVBasis::evaluateOperatorDense(const ScalarGSQOneElectronOperator<double>&): The number of orbitals of this ONV basis and the given one-electron operator are incompatible.");
        }

        // Initialize a container for the dense matrix representation, and fill it with the general evaluation function.
        MatrixRepresentationEvaluationContainer<SquareMatrix<Scalar>> container {this->dimension()};
        this->evaluate<SquareMatrix<Scalar>>(f, container);

        return container.evaluation();
    }


    /**
     *  Calculate the dense matrix representation of a generalized Hamiltonian in this ONV basis.
     *
     *  @tparam Scalar          The scalar representation of a one-electron parameter: real or complex.
     *
     *  @param hamiltonian      A generalized Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the Hamiltonian.
     */
    template <typename Scalar>
    SquareMatrix<Scalar> evaluateOperatorDense(const GSQHamiltonian<Scalar>& hamiltonian) const {

        if (hamiltonian.numberOfOrbitals() != this->numberOfOrbitals()) {
            throw std::invalid_argument("SpinUnresolvedSelectedONVBasis::evaluateOperatorDense(const GSQHamiltonian<double>&): The number of orbitals of this ONV basis and the given Hamiltonian are incompatible.");
        }

        // Initialize a container for the dense matrix representation, and fill it with the general evaluation function.
        MatrixRepresentationEvaluationContainer<SquareMatrix<Scalar>> container {this->dimension()};
        this->evaluate<SquareMatrix<Scalar>>(hamiltonian, container);

        return container.evaluation();
    }


    /*
     *  MARK: Operator evaluations - general implementations - containers
     */

    /**
     *  Calculate the matrix representation of a generalized one-electron operator in this ONV basis and emplace it in the given container.
     *
     *  @tparam Matrix                      The type of matrix used to store the evaluations.
     *  @tparam Scalar                      The scalar representation of a one-electron parameter: real or complex.
     *
     *  @param f_op                         A generalized one-electron operator expressed in an orthonormal spinor basis.
     *  @param container                    A specialized container for emplacing evaluations/matrix elements.
     */
    template <typename Matrix, typename Scalar>
    void evaluate(const ScalarGSQOneElectronOperator<Scalar>& f_op, MatrixRepresentationEvaluationContainer<Matrix>& container) const {

        // Prepare some variables.
        const auto& f = f_op.parameters();
        const auto dim = this->dimension();


        // Loop over all bra indices I.
        for (; !container.isFinished(); container.increment()) {
            auto onv_I = this->onvWithIndex(container.index);

            // Calculate the diagonal elements (I = J).
            for (size_t p = 0; p < this->numberOfOrbitals(); p++) {
                if (onv_I.isOccupied(p)) {
                    container.addRowwise(container.index, f(p, p));  // This emplaces F(I,I).
                }
            }


            // Calculate the off-diagonal elements (I != J), by going over all other ket ONVs J.
            for (size_t J = container.index + 1; J < dim; J++) {
                auto onv_J = this->onvWithIndex(J);

                // If I and J are only 1 excitation away, they can couple through the operator.
                if (onv_I.countNumberOfExcitations(onv_J) == 1) {
                    auto p = onv_I.findDifferentOccupations(onv_J)[0];  // The orbital that is occupied in I, but not in J.
                    auto q = onv_J.findDifferentOccupations(onv_I)[0];  // The orbital that is occupied in J, but not in I.
                    // We're sure that there is only 1 element in the vectors above.

                    // Calculate the total sign and emplace the correct value in the container.
                    const auto sign = static_cast<double>(onv_I.operatorPhaseFactor(p) * onv_J.operatorPhaseFactor(q));
                    const auto value = sign * f(p, q);

                    container.addColumnwise(J, value);           // This emplaces F(I,J).
                    container.addRowwise(J, GQCP::conj(value));  // This emplaces F(J,I).
                }
            }
        }
    }


    /**
     *  Calculate the matrix representation of a generalized Hamiltonian in this ONV basis and emplace it in the given container.
     *
     *  @tparam Matrix                      The type of matrix used to store the evaluations.
     *  @tparam Scalar                      The scalar representation of a Hamiltonian element: real or complex.
     *
     *  @param hamiltonian                  A generalized Hamiltonian expressed in an orthonormal spinor basis.
     *  @param container                    A specialized container for emplacing evaluations/matrix elements.
     */
    template <typename Matrix, typename Scalar>
    void evaluate(const GSQHamiltonian<Scalar>& hamiltonian, MatrixRepresentationEvaluationContainer<Matrix>& container) const {

        // Prepare some variables.
        const size_t dim = this->dimension();
        const size_t M = this->numberOfOrbitals();

        const auto& h = hamiltonian.core().parameters();
        const auto& g = hamiltonian.twoElectron().parameters();


        // Loop over all bra indices I.
        for (; !container.isFinished(); container.increment()) {
            auto onv_I = this->onvWithIndex(container.index);
            auto& occupied_indices_I = onv_I.occupiedIndices();

            // Calculate the diagonal elements (I = J).
            for (const auto p : occupied_indices_I) {
                container.addRowwise(container.index, h(p, p));  // This emplaces F(I,I).

                for (const auto q : occupied_indices_I) {
                    if (p != q) {
                        container.addRowwise(container.index, 0.5 * (g(p, p, q, q) - g(p, q, q, p)));
                    }
                }
            }


            // Calculate the off-diagonal elements (I != J), by going over all other ket ONVs J.
            for (size_t J = container.index + 1; J < dim; J++) {
                auto onv_J = this->onvWithIndex(J);

                // If I and J are only 1 excitation away, they can couple through the Hamiltonian through both the one- and two-electron parts.
                if (onv_I.countNumberOfExcitations(onv_J) == 1) {

                    // The one-electron part.
                    auto w = onv_I.findDifferentOccupations(onv_J)[0];  // The orbital that is occupied in I, but not in J.
                    auto x = onv_J.findDifferentOccupations(onv_I)[0];  // The orbital that is occupied in J, but not in I.
                    // We're sure that there is only 1 element in the vectors above.

                    // Calculate the total sign and emplace the correct value in the container.
                    const auto sign = static_cast<double>(onv_I.operatorPhaseFactor(w) * onv_J.operatorPhaseFactor(x));
                    const auto value = sign * h(w, x);

                    container.addColumnwise(J, value);           // This emplaces F(I,J).
                    container.addRowwise(J, GQCP::conj(value));  // This emplaces F(J,I).


                    // The two-electron part.
                    for (size_t p = 0; p < M; p++) {
                        if (onv_I.isOccupied(p) && onv_J.isOccupied(p)) {  // `p` must be occupied in I and J.
                            if ((p != w) && (p != x)) {                    // We can't annihilate twice on the indices `w` or `x`.

                                const auto value = 0.5 * sign * (g(w, x, p, p) - g(w, p, p, x) + g(p, p, w, x) - g(p, x, w, p));

                                container.addColumnwise(J, value);           // This emplaces G(I,J).
                                container.addRowwise(J, GQCP::conj(value));  // This emplaces G(J,I).
                            }
                        }
                    }
                }

                // If I and J are 2 excitations away, they can couple through the Hamiltonian through the two-electron part.
                else if (onv_I.countNumberOfExcitations(onv_J) == 2) {

                    auto occupied_in_I = onv_I.findDifferentOccupations(onv_J);  // The orbitals that are occupied in J, but not in J.
                    const auto w = occupied_in_I[0];
                    const auto x = occupied_in_I[1];

                    auto occupied_in_J = onv_J.findDifferentOccupations(onv_I);  // The orbitals that are occupied in I, but not in J.
                    const auto y = occupied_in_J[0];
                    const auto z = occupied_in_J[1];

                    const auto sign = static_cast<double>(onv_I.operatorPhaseFactor(w) * onv_I.operatorPhaseFactor(x) * onv_J.operatorPhaseFactor(y) * onv_J.operatorPhaseFactor(z));
                    const auto value = 0.5 * sign * (g(x, z, w, y) - g(w, z, x, y) + g(w, y, x, z) - g(x, y, w, z));

                    container.addColumnwise(J, value);           // This emplaces G(I,J).
                    container.addRowwise(J, GQCP::conj(value));  // This emplaces G(J,I).
                }
            }  // Loop over ket addresses J > I.
        }      // Loop over bra addresses I.
    }
};


}  // namespace GQCP
