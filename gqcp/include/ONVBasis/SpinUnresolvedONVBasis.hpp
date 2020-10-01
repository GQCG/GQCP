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
#include "ONVBasis/BaseONVBasis.hpp"
#include "ONVBasis/ONVManipulator.hpp"
#include "ONVBasis/ONVPath.hpp"

#include <functional>


namespace GQCP {


/**
 *  The full spin-unresolved ONV basis for a number of spinors and number of electrons.
 */
class SpinUnresolvedONVBasis:
    public BaseONVBasis,
    public ONVManipulator<SpinUnresolvedONVBasis> {
private:
    std::vector<std::vector<size_t>> vertex_weights;  // vertex_weights of the addressing scheme

public:
    using ONV = SpinUnresolvedONV;


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  The default constructor
     */
    SpinUnresolvedONVBasis() = default;

    /**
     *  @param M            the number of spinors
     *  @param N            the number of electrons
     */
    SpinUnresolvedONVBasis(const size_t M, const size_t N);


    // DESTRUCTOR

    /**
     *  The default destructor
     */
    ~SpinUnresolvedONVBasis() override = default;


    // STATIC PUBLIC METHODS

    /**
     *  @param M            the number of spinors
     *  @param N            the number of electrons
     *
     *  @return the dimension of a spin-unresolved ONV basis with the given number of spinors and electrons
     */
    static size_t calculateDimension(const size_t M, const size_t N);


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param representation      a representation of an spin-unresolved ONV
     *
     *  @return the address (i.e. the ordering number) of the given spin-unresolved ONV
     */
    size_t addressOf(const size_t representation) const override;

    /**
     *  @param onv       the spin-unresolved ONV
     *
     *  @return the number of ONVs (with a larger address) this spin-unresolved ONV would couple with given a one electron operator
     */
    size_t countOneElectronCouplings(const SpinUnresolvedONV& onv) const override;

    /**
     *  @return the number of non-zero (non-diagonal) couplings of a one electron coupling scheme in the spin-unresolved ONV basis
     */
    size_t countTotalOneElectronCouplings() const override;

    /**
     *  @return the number of non-zero (non-diagonal) couplings of a two electron coupling scheme in the spin-unresolved ONV basis
     */
    size_t countTotalTwoElectronCouplings() const override;

    /**
     *  @param onv       the spin-unresolved ONV
     *
     *  @return the number of ONVs (with a larger address) this spin-unresolved ONV would couple with given a two electron operator
     */
    size_t countTwoElectronCouplings(const SpinUnresolvedONV& onv) const override;

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the spin-unresolved ONV basis
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarSQOneElectronOperator<double>& one_op, const bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the spin-unresolved ONV basis
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarSQTwoElectronOperator<double>& two_op, const bool diagonal_values) const override;

    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values              bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the spin-unresolved ONV basis
     */
    SquareMatrix<double> evaluateOperatorDense(const SQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const override;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the spin-unresolved ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const override;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the spin-unresolved ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarSQTwoElectronOperator<double>& two_op) const override;

    /**
     *  Evaluate the diagonal of the Hamiltonian
     *
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the spin-unresolved ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the spin-unresolved ONV basis
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQOneElectronOperator<double>& one_op, const bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the spin-unresolved ONV basis
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQTwoElectronOperator<double>& two_op, const bool diagonal_values) const override;

    /**
     *  Evaluate the Hamiltonian in a sparse matrix
     *
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values              bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the spin-unresolved ONV basis
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const SQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const override;

    /**
     *  @param representation       a representation of an spin-unresolved ONV
     *
     *  @return the next bitstring permutation in the spin-unresolved ONV basis
     *
     *      Examples:
     *          011 -> 101
     *          101 -> 110
     */
    size_t nextPermutationOf(const size_t representation) const override;

    /**
      *  Calculate unsigned representation for a given address
      *
      *  @param address                 the address of the representation is calculated
      *
      *  @return unsigned representation of the address
      */
    size_t representationOf(const size_t address) const override;


    // PUBLIC METHODS

    using ONVManipulator<SpinUnresolvedONVBasis>::addressOf;

    /**
     *  @param p       The orbital index p
     *  @param n       The electron index n
     *
     *  @return the arc weight of the arc starting in the given vertex (p, n) (reference: Helgaker 11.8.6)
     */
    size_t arcWeight(const size_t p, const size_t n) const { return this->vertexWeight(p, n + 1); }

    /**
     *  Calculates sigma(pq) + sigma(qp)'s: all one-electron couplings for each annihilation-creation pair in the (spin) spin-unresolved ONV basis
     *  and stores them in sparse matrices for each pair combination
     *
     *  @return vector of sparse matrices containing the one-electron couplings for the (spin) spin-unresolved ONV basis
     *      Ordered as: sigma(00), sigma(01) + sigma(10), sigma(02)+ sigma(20), ...
     */
    std::vector<Eigen::SparseMatrix<double>> calculateOneElectronCouplings() const;

    /**
     *  Evaluate the operator in a given evaluation iterator in the spin-unresolved ONV basis
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param one_op                        the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param evaluation_iterator           evaluation iterator to which the evaluations are added
     *  @param diagonal_values               bool to indicate if diagonal values will be calculated
     */
    template <typename _Matrix>
    void evaluate(const ScalarSQOneElectronOperator<double>& one_op, MatrixRepresentationEvaluationContainer<_Matrix>& evaluation_iterator, const bool diagonal_values) const {

        const auto& one_op_par = one_op.parameters();

        const size_t K = this->numberOfOrbitals();
        const size_t N = this->numberOfElectrons();
        const size_t dim = this->dimension();

        SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0

        for (; !evaluation_iterator.isFinished(); evaluation_iterator.increment()) {  // I loops over all the addresses of the onv
            for (size_t e1 = 0; e1 < N; e1++) {                                       // e1 (electron 1) loops over the (number of) electrons

                size_t p = onv.occupationIndexOf(e1);  // retrieve the index of a given electron
                // remove the weight from the initial address I, because we annihilate
                size_t address = evaluation_iterator.index - this->vertexWeight(p, e1 + 1);

                if (diagonal_values) {
                    evaluation_iterator.addRowwise(evaluation_iterator.index, one_op_par(p, p));
                }

                // The e2 iteration counts the number of encountered electrons for the creation operator
                // We only consider greater addresses than the initial one (because of symmetry)
                // Hence we only count electron after the annihilated electron (e1)
                size_t e2 = e1 + 1;
                size_t q = p + 1;

                int sign_e2 = 1;
                // perform a shift
                this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);

                while (q < K) {
                    size_t J = address + this->vertexWeight(q, e2);
                    double value = sign_e2 * one_op_par(p, q);
                    evaluation_iterator.addColumnwise(J, value);
                    evaluation_iterator.addRowwise(J, value);

                    q++;  // go to the next orbital

                    // perform a shift
                    this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
                }  //  (creation)
            }      // e1 loop (annihilation)
            // Prevent last permutation
            if (evaluation_iterator.index < dim - 1) {
                this->transformONVToNextPermutation(onv);
            }
        }
    }

    template <typename Representation>
    Representation evaluate_new(const ScalarSQOneElectronOperator<double>& one_op) const {

        const auto& one_op_par = one_op.parameters();
        const size_t dim = this->dimension();

        MatrixRepresentationEvaluationContainer<Representation> ONV_iterator {dim};
        SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0

        for (; !ONV_iterator.isFinished(); ONV_iterator.increment()) {  // loops over all possible ONVs

            for (size_t e1 = 0; e1 < N; e1++) {  // loop over electrons that can be annihilated

                ONVPath<SpinUnresolvedONVBasis> onv_path {*this, onv};

                std::cout << "\n\n"
                          << onv_path.address() << "    " << onv_path.orbitalIndex() << "    " << onv_path.electronIndex() << std::endl;

                size_t q = onv.occupationIndexOf(e1);  // retrieve orbital index of the electron

                onv_path.annihilate(q, e1);
                std::cout << "\n"
                          << "a: " << onv_path.address() << "    " << onv_path.orbitalIndex() << "    " << onv_path.electronIndex() << std::endl;

                while (!onv_path.isFinished()) {

                    onv_path.leftTranslateUntilVertical();
                    std::cout << "LTUV: " << onv_path.address() << "    " << onv_path.orbitalIndex() << "    " << onv_path.electronIndex() << std::endl;
                    onv_path.create();
                    std::cout << "C: " << onv_path.address() << "    " << onv_path.orbitalIndex() << "    " << onv_path.electronIndex() << std::endl;

                    const auto address = onv_path.address();

                    const auto h_pq = onv_path.sign() * one_op_par(onv_path.orbitalIndex() - 1, q);

                    ONV_iterator.addColumnwise(address, h_pq);
                    ONV_iterator.addRowwise(address, h_pq);

                    if () {
                        onv_path.annihilate(onv_path.orbitalIndex() - 1, onv_path.electronIndex() - 1);
                        std::cout << "A: " << onv_path.address() << "    " << onv_path.orbitalIndex() << "    " << onv_path.electronIndex() << std::endl;
                    }
                }
            }
            // Prevent last ONV since there is no possibility for an electron to be annihilated anymore.
            if (ONV_iterator.index < dim - 1) {
                this->transformONVToNextPermutation(onv);
            }
        }

        return ONV_iterator.evaluation();
    }


    /**
     *  Evaluate a one-electron operator in this spin-unresolved ONV basis.
     * 
     *  @tparam Representation              The matrix representation that is used for storing the result. Essentially, any type that can be used in MatrixRepresentationEvaluationContainer<Representation>.
     * 
     *  @param one_op                       A one-electron operator in an orthonormal orbital basis.
     *  @param should_calculate_diagonal    If diagonal values should be calculated.
     * 
     *  @return The matrix representation of the given one-electron operator.
     */
    template <typename Representation>
    Representation evaluate(const ScalarSQOneElectronOperator<double>& one_op, const bool diagonal_values) const {

        const auto& one_op_par = one_op.parameters();

        MatrixRepresentationEvaluationContainer<Representation> evaluation_iterator {this->dimension()};

        const size_t K = this->numberOfOrbitals();
        const size_t N = this->numberOfElectrons();
        const size_t dim = this->dimension();

        SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0

        for (; !evaluation_iterator.isFinished(); evaluation_iterator.increment()) {  // I loops over all the addresses of the onv
            for (size_t e1 = 0; e1 < N; e1++) {                                       // e1 (electron 1) loops over the (number of) electrons

                size_t p = onv.occupationIndexOf(e1);  // retrieve the index of a given electron
                // remove the weight from the initial address I, because we annihilate
                size_t address = evaluation_iterator.index - this->vertexWeight(p, e1 + 1);

                if (diagonal_values) {
                    evaluation_iterator.addRowwise(evaluation_iterator.index, one_op_par(p, p));
                }

                // The e2 iteration counts the number of encountered electrons for the creation operator
                // We only consider greater addresses than the initial one (because of symmetry)
                // Hence we only count electron after the annihilated electron (e1)
                size_t e2 = e1 + 1;
                size_t q = p + 1;

                int sign_e2 = 1;
                // perform a shift
                this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);

                while (q < K) {
                    size_t J = address + this->vertexWeight(q, e2);
                    double value = sign_e2 * one_op_par(p, q);
                    evaluation_iterator.addColumnwise(J, value);
                    evaluation_iterator.addRowwise(J, value);

                    q++;  // go to the next orbital

                    // perform a shift
                    this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
                }  //  (creation)
            }      // e1 loop (annihilation)
            // Prevent last permutation
            if (evaluation_iterator.index < dim - 1) {
                this->transformONVToNextPermutation(onv);
            }
        }
        return evaluation_iterator.evaluation();
    }


    /**
     *  Evaluate the operator in a given evaluation iterator in the spin-unresolved ONV basis
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param two_op                        the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param evaluation_iterator           evaluation iterator to which the evaluations are added
     *  @param diagonal_values               bool to indicate if diagonal values will be calculated
     */
    template <typename _Matrix>
    void evaluate(const ScalarSQTwoElectronOperator<double>& two_op, MatrixRepresentationEvaluationContainer<_Matrix>& evaluation_iterator, const bool diagonal_values) const {

        // Calling this combined method for both the one- and two-electron operator does not affect the performance, hence we avoid writing more code by plugging a zero operator in the combined method
        evaluate(ScalarSQOneElectronOperator<double> {this->M}, two_op, evaluation_iterator, diagonal_values);
    }


    /**
     *  Evaluate the operators in a given evaluation iterator in the spin-unresolved ONV basis
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param one_op                        the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param two_op                        the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param evaluation_iterator           evaluation iterator to which the evaluations are added
     *  @param diagonal_values               bool to indicate if diagonal values will be calculated
     */
    template <typename _Matrix>
    void evaluate(const ScalarSQOneElectronOperator<double>& one_op, const ScalarSQTwoElectronOperator<double>& two_op, MatrixRepresentationEvaluationContainer<_Matrix>& evaluation_iterator, const bool diagonal_values) const {

        const auto& two_op_par = two_op.parameters();

        const size_t K = this->numberOfOrbitals();
        const size_t N = this->numberOfElectrons();
        const size_t dim = this->dimension();

        ScalarSQOneElectronOperator<double> k = two_op.effectiveOneElectronPartition() + one_op;
        const auto& k_par = k.parameters();

        SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0

        for (; !evaluation_iterator.isFinished(); evaluation_iterator.increment()) {  // I loops over all addresses in the spin-unresolved ONV basis
            if (evaluation_iterator.index > 0) {
                this->transformONVToNextPermutation(onv);
            }
            int sign1 = -1;                      // start with -1 because we flip at the start of the annihilation (so we start at 1, followed by:  -1, 1, ...)
            for (size_t e1 = 0; e1 < N; e1++) {  // A1 (annihilation 1)

                sign1 *= -1;
                size_t p = onv.occupationIndexOf(e1);  // retrieve the index of a given electron
                size_t address = evaluation_iterator.index - this->vertexWeight(p, e1 + 1);

                // Strictly diagonal values
                if (diagonal_values) {
                    evaluation_iterator.addRowwise(evaluation_iterator.index, k_par(p, p));
                    for (size_t q = 0; q < K; q++) {  // q loops over SOs
                        if (onv.isOccupied(q)) {
                            evaluation_iterator.addRowwise(evaluation_iterator.index, 0.5 * two_op_par(p, p, q, q));
                        } else {
                            evaluation_iterator.addRowwise(evaluation_iterator.index, 0.5 * two_op_par(p, q, q, p));
                        }
                    }
                }

                size_t address1 = address;
                size_t e2 = e1;
                size_t q = p;

                /**
                 *  A1 > C1 (annihlation 1 > creation 1)
                 */
                int sign2 = sign1;
                q--;
                e2--;
                this->shiftUntilPreviousUnoccupiedOrbital<1>(onv, address1, q, e2, sign2);
                while (q != -1) {

                    size_t address2 = address1 + this->vertexWeight(q, e2 + 2);

                    /**
                     *  C2 > A2
                     */
                    int sign3 = sign1;
                    for (size_t e3 = e1 + 1; e3 < N; e3++) {
                        sign3 *= -1;  // initial sign3 = sign of the annihilation, with one extra electron(from crea) = *-1
                        size_t r = onv.occupationIndexOf(e3);
                        size_t address3 = address2 - this->vertexWeight(r, e3 + 1);

                        size_t e4 = e3 + 1;
                        size_t s = r + 1;

                        int sign4 = sign3;
                        this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                        while (s < K) {
                            size_t J = address3 + this->vertexWeight(s, e4);
                            int signev = sign1 * sign2 * sign3 * sign4;
                            double value = signev * 0.5 * (two_op_par(p, q, r, s) + two_op_par(r, s, p, q) - two_op_par(p, s, r, q) - two_op_par(r, q, p, s));


                            evaluation_iterator.addColumnwise(J, value);
                            evaluation_iterator.addRowwise(J, value);

                            s++;
                            this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);
                        }
                    }
                    q--;
                    this->shiftUntilPreviousUnoccupiedOrbital<1>(onv, address1, q, e2, sign2);
                }

                e2 = e1 + 1;
                q = p + 1;
                sign2 = sign1;
                this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign2);

                /**
                 *  A1 < C1
                 */
                while (q < K) {

                    address1 = address + this->vertexWeight(q, e2);

                    /**
                     *  A2 > C1
                     */
                    int sign3 = sign2;
                    for (size_t e3 = e2; e3 < N; e3++) {
                        sign3 *= -1;  // -1 cause we created electron (creation) sign of A is now the that of C *-1
                        size_t r = onv.occupationIndexOf(e3);
                        size_t address3 = address1 - this->vertexWeight(r, e3 + 1);

                        size_t e4 = e3 + 1;
                        size_t s = r + 1;
                        int sign4 = sign3;
                        this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                        while (s < K) {
                            size_t J = address3 + this->vertexWeight(s, e4);
                            int signev = sign1 * sign2 * sign3 * sign4;

                            double value = signev * 0.5 * (two_op_par(p, q, r, s) + two_op_par(r, s, p, q) - two_op_par(r, q, p, s) - two_op_par(p, s, r, q));

                            evaluation_iterator.addColumnwise(J, value);
                            evaluation_iterator.addRowwise(J, value);

                            s++;  // go to the next orbital
                            this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                        }  // (creation)
                    }

                    size_t r = q;
                    sign3 = sign2;
                    size_t address1c = address1;

                    /**
                     *  A2 < C1, (A2 > A1)
                     */
                    for (size_t e3 = e2 - 1; e3 > e1; e3--) {
                        sign3 *= -1;
                        size_t e4 = e2;
                        address1c += this->vertexWeight(r, e3) - this->vertexWeight(r, e3 + 1);
                        r = onv.occupationIndexOf(e3);
                        size_t address2 = address1c - this->vertexWeight(r, e3);
                        int sign4 = sign2;
                        size_t s = q + 1;
                        this->shiftUntilNextUnoccupiedOrbital<1>(onv, address2, s, e4, sign4);
                        while (s < K) {

                            size_t J = address2 + this->vertexWeight(s, e4);

                            int signev = sign1 * sign2 * sign3 * sign4;
                            double value = signev * 0.5 * (two_op_par(p, q, r, s) + two_op_par(r, s, p, q) - two_op_par(r, q, p, s) - two_op_par(p, s, r, q));

                            evaluation_iterator.addColumnwise(J, value);
                            evaluation_iterator.addRowwise(J, value);

                            s++;
                            this->shiftUntilNextUnoccupiedOrbital<1>(onv, address2, s, e4, sign4);
                        }
                    }

                    /**
                     *  A2 = C1
                     */
                    int signev = sign2 * sign1;

                    double value_I = k_par(p, q);  // cover the one electron calculations

                    for (size_t s = 0; s < K; s++) {
                        if (!onv.isOccupied(s)) {
                            value_I += 0.5 * (two_op_par(p, s, s, q));
                        } else {
                            value_I += 0.5 * (two_op_par(s, s, p, q) - two_op_par(s, q, p, s) + two_op_par(p, q, s, s));
                        }
                    }

                    value_I *= signev;

                    q++;

                    evaluation_iterator.addColumnwise(address1, value_I);
                    evaluation_iterator.addRowwise(address1, value_I);

                    this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign2);
                }
            }
        }
    }


    /**
     *  Evaluate the matrix-vector product of a one-electron operator
     *
     *  @param one_op                       the one-electron operator expressed in an orthonormal basis
     *  @param x                            the vector of the matrix-vector product
     *  @param diagonal                     the diagonal of the matrix representation of the operator inside the spin-unresolved ONV basis
     *
     *  @return a vector that is equal to the matrix-vector product of the one-electron operator's matrix representation and the given vector
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Evaluate a two electron operator in a matrix vector product
     *
     *  @param two_op                       the two electron operator expressed in an orthonormal basis
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in the spin-unresolved ONV basis
     *
     *  @return a vector that is equal to the matrix-vector product of the two-electron operator's matrix representation and the given vector
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarSQTwoElectronOperator<double>& two_op, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Evaluate the Hamiltonian in a matrix vector product
     *
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in the spin-unresolved ONV basis
     *
     *  @return a vector that is equal to the matrix-vector product of the Hamiltonian's matrix representation and the given vector
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Iterate over all ONVs in this ONV basis and apply the given callback function.
     * 
     *  @param callback             the function to be applied in every iteration. Its arguments are a spin-unresolved ONV and its corresponding addresses
     */
    void forEach(const std::function<void(const SpinUnresolvedONV&, const size_t)>& callback) const;

    /**
     *  Find the next unoccupied orbital in a given spin-unresolved ONV,
     *  update the electron count, orbital index,
     *  and update the address by calculating a shift
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for previously annihilated electrons
     *
     *  @tparam T        the number of previously annihilated electrons
     *
     *  @param address   the address which is updated
     *  @param onv       the spin-unresolved ONV for which we search the next unnocupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     */
    template <int T>
    void shiftUntilNextUnoccupiedOrbital(const SpinUnresolvedONV& onv, size_t& address, size_t& q, size_t& e) const {

        // Test whether the current orbital index is occupied
        while (e < this->N && q == onv.occupationIndexOf(e)) {

            // Take the difference of vertex weights for the encountered electron weights to that of a vertex weight path with "a" fewer electrons
            // +1 is added to the electron index, because of how the addressing scheme is arrayed.
            address += this->vertexWeight(q, e + 1 - T) - this->vertexWeight(q, e + 1);

            // move to the next electron and orbital
            e++;
            q++;
        }
    }


    /**
     *  Find the next unoccupied orbital in a given spin-unresolved ONV,
     *  update the electron count, orbital index, sign,
     *  and update the address by calculating a shift
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for previously annihilated electrons
     *
     *  @tparam T        the number of previously annihilated electrons
     *
     *  @param address   the address which is updated
     *  @param onv       the spin-unresolved ONV for which we search the next unnocupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     *  @param sign      the sign which is flipped for each iteration
     */
    template <int T>
    void shiftUntilNextUnoccupiedOrbital(const SpinUnresolvedONV& onv, size_t& address, size_t& q, size_t& e, int& sign) const {

        // Test whether the current orbital index is occupied
        while (e < this->N && q == onv.occupationIndexOf(e)) {

            // Take the difference of vertex weights for the encountered electron weights to that of a vertex weight path with "a" fewer electrons
            // +1 is added to the electron index, because of how the addressing scheme is arrayed.
            address += this->vertexWeight(q, e + 1 - T) - this->vertexWeight(q, e + 1);

            // move to the next electron and orbital
            e++;
            q++;
            sign *= -1;
        }
    }


    /**
     *  Find the previous unoccupied orbital in a given spin-unresolved ONV,
     *  update the electron count, orbital index, sign,
     *  and update the address by calculating a shift
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for newly created electrons
     *
     *  @tparam T        the number of newly created electrons
     *
     *  @param address   the address which is updated
     *  @param onv       the spin-unresolved ONV for which we search the next unoccupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     *  @param sign      the sign which is flipped for each iteration
     */
    template <int T>
    void shiftUntilPreviousUnoccupiedOrbital(const SpinUnresolvedONV& onv, size_t& address, size_t& q, size_t& e, int& sign) const {

        // Test whether the current orbital index is occupied
        while (e != -1 && q == onv.occupationIndexOf(e)) {

            int shift = static_cast<int>(this->vertexWeight(q, e + 1 + T)) - static_cast<int>(this->vertexWeight(q, e + 1));
            address += shift;

            e--;
            q--;
            sign *= -1;
        }
    }


    /**
     *  @param p            the orbital index
     *  @param n            the electron index
     * 
     *  @return the vertex weight for the given indices
     */
    size_t vertexWeight(const size_t p, const size_t n) const { return this->vertex_weights[p][n]; }

    /**
     *  @return all the vertex weights for this ONV basis
     */
    const std::vector<std::vector<size_t>>& vertexWeights() const { return this->vertex_weights; }
};


}  // namespace GQCP
