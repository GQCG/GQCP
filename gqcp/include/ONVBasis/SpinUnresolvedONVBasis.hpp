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
#include "ONVBasis/ONVPath.hpp"
#include "ONVBasis/SpinUnresolvedONV.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/GSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/PureUSQTwoElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperatorComponent.hpp"
#include "Utilities/complex.hpp"

#include <functional>


namespace GQCP {


/**
 *  The full spin-unresolved ONV basis for a number of spinors/spin-orbitals and number of electrons.
 */
class SpinUnresolvedONVBasis {
private:
    // The number of spinors/spin-orbitals.
    size_t M;

    // The number of electrons, i.e. the number of occupied spinors/spin-orbitals.
    size_t N;

    // The vertex weights corresponding to the addressing scheme for a full spin-unresolved ONV basis. This addressing scheme is taken from Helgaker, Jørgensen, Olsen (2000).
    std::vector<std::vector<size_t>> vertex_weights;

public:
    // The ONV that is naturally related to a full spin-unresolved ONV basis. See also `ONVPath`.
    using ONV = SpinUnresolvedONV;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param M            The number of spinors/spin-orbitals.
     *  @param N            The number of electrons, i.e. the number of occupied spinors/spin-orbitals.
     */
    SpinUnresolvedONVBasis(const size_t M, const size_t N);

    /**
     *  The default constructor.
     */
    SpinUnresolvedONVBasis() = default;


    /*
     *  MARK: Basic information
     */

    /**
     *  Calculate the dimension of the full spin-unresolved ONV basis with a given number of spinors/spin-orbitals and number of electrons.
     * 
     *  @param M            The number of spinors/spin-orbitals.
     *  @param N            The number of electrons, i.e. the number of occupied spinors/spin-orbitals.
     *
     *  @return The dimension of a spin-unresolved ONV basis.
     */
    static size_t calculateDimension(const size_t M, const size_t N);

    /**
     *  @return The dimension of this ONV basis.
     */
    size_t dimension() const { return SpinUnresolvedONVBasis::calculateDimension(this->numberOfOrbitals(), this->numberOfElectrons()); }

    /**
     *  @return The number of electrons, i.e. the number of occupied spinors/spin-orbitals.
     */
    size_t numberOfElectrons() const { return this->N; }

    /**
     *  @return The number of spinors/spin-orbitals.
     */
    size_t numberOfOrbitals() const { return this->M; }


    /**
     *  MARK: Addressing scheme, address calculations and ONV manipulations
     */

    /**
     *  Access the arc weight of an arc in the addressing scheme of this ONV basis. The addressing scheme is taken from Helgaker, Jørgensen, Olsen (2000).
     * 
     *  @param p            The orbital index.
     *  @param n            The electron index.
     *
     *  @return The arc weight of the arc starting at the given vertex (p, n).
     */
    size_t arcWeight(const size_t p, const size_t n) const;

    /**
     *  @param p            The orbital index.
     *  @param n            The electron index.
     * 
     *  @return The vertex weight related to the given indices (p,n).
     */
    size_t vertexWeight(const size_t p, const size_t n) const { return this->vertex_weights[p][n]; }

    /**
     *  @return All the vertex weights for this ONV basis, stored as a vector of vectors. The outer axis represents the orbital indices, the inner axis represents the electron indices.
     */
    const std::vector<std::vector<size_t>>& vertexWeights() const { return this->vertex_weights; }

    /**
     *  Calculate the address (i.e. the ordering number) of an unsigned representation of a spin-unresolved ONV.
     * 
     *  @param representation      The unsigned representation of a spin-unresolved ONV.
     *
     *  @return The address corresponding to the unsigned representation of a spin-unresolved ONV.
     */
    size_t addressOf(const size_t representation) const;

    /**
     *  Calculate the address (i.e. the ordering number) of a spin-unresolved ONV.
     * 
     *  @param onv          The spin-unresolved ONV.
     *
     *  @return The address (i.e. the ordering number) of the given spin-unresolved ONV.
     */
    size_t addressOf(const SpinUnresolvedONV& onv) const { return this->addressOf(onv.unsignedRepresentation()); }

    /**
     *  Calculate the next allowed unsigned representation of a spin-unresolved ONV in this ONV basis.
     * 
     *  @param representation       A representation of a spin-unresolved ONV.
     *
     *  @return The next allowed unsigned representation of a spin-unresolved ONV in this ONV basis.
     *
     *  @example
     *          011 -> 101
     *          101 -> 110
     */
    size_t nextPermutationOf(const size_t representation) const;

    /**
     *  Calculate the unsigned representation of a spin-unresolved ONV that corresponds to the address/ordering number in this ONV basis.
     *
     *  @param address                 The address/ordering number of a spin-unresolved ONV in this ONV basis.
     *
     *  @return The unsigned representation of a spin-unresolved ONV that corresponds to the address/ordering number in this ONV basis.
     */
    size_t representationOf(const size_t address) const;

    /**
     *  Create the ONV that corresponds to the given address in this ONV basis.
     * 
     *  @param address                 The address/ordering number of a spin-unresolved ONV in this ONV basis.
     *
     *  @return The ONV that corresponds to the given address in this ONV basis.
     */
    SpinUnresolvedONV constructONVFromAddress(const size_t address) const;

    /**
     *  Modify a `SpinResolvedONV` to the next allowed ONV in this ONV basis.
     * 
     *  @param onv      A spin-resolved ONV.
     */
    void transformONVToNextPermutation(SpinUnresolvedONV& onv) const;

    /**
     *  Modify a `SpinResolvedONV` to the one with the given address in this ONV basis.
     *
     *  @param onv          A spin-resolved ONV.
     *  @param address      The target address in this ONV basis.
     */
    void transformONVCorrespondingToAddress(SpinUnresolvedONV& onv, const size_t address) const;


    /*
     *  MARK: Couplings
     */

    /**
     *  Calculate the number of ONVs (with a larger address) that a given spin-unresolved ONV would couple with in this ONV basis, through a one-electron operator.
     * 
     *  @param onv          The spin-unresolved ONV.
     *
     *  @return The number of ONVs (with a larger address) that a given spin-unresolved ONV would couple with in this ONV basis, through a one-electron operator.
     */
    size_t countOneElectronCouplings(const SpinUnresolvedONV& onv) const;

    /**
     *  @return The total number of non-zero and non-diagonal couplings of a one-electron operator in this ONV basis.
     */
    size_t countTotalOneElectronCouplings() const;

    /**
     *  @return The total number of non-zero and non-diagonal couplings of a two-electron operator in this ONV basis.
     */
    size_t countTotalTwoElectronCouplings() const;

    /**
     *  Calculate the number of ONVs (with a larger address) that a given spin-unresolved ONV would couple with in this ONV basis, through a two-electron operator.
     * 
     *  @param onv          The spin-unresolved ONV.
     *
     *  @return The number of ONVs (with a larger address) that a given spin-unresolved ONV would couple with in this ONV basis, through a two-electron operator.
     */
    size_t countTwoElectronCouplings(const SpinUnresolvedONV& onv) const;

    /**
     *  Calculate all one-electron coupling elements for this spin-unresolved ONV basis. These are all the intermediate matrices sigma(pq)_{IJ}, as defined by Helgaker, Jørgensen, Olsen (2000).
     *
     *  @return A vector of sparse matrices containing the one-electron coupling elements for this spin-unresolved ONV basis. The elements of this vector are ordered through the one-electron excitation (pq) in ascending order: sigma(00), sigma(01) + sigma(10), sigma(02)+ sigma(20), ...
     */
    std::vector<Eigen::SparseMatrix<double>> calculateOneElectronCouplings() const;


    /**
     *  MARK: Iterating
     */

    /**
     *  Iterate over all ONVs in this ONV basis and apply the given callback function.
     * 
     *  @param callback            The function to be applied in every iteration. Its supplied arguments are a spin-unresolved ONV and its corresponding addresses.
     */
    void forEach(const std::function<void(const SpinUnresolvedONV&, const size_t)>& callback) const;


    /*
     *  MARK: Dense generalized operator evaluations
     */

    /**
     *  Calculate the dense matrix representation of a generalized one-electron operator in this ONV basis.
     *
     *  @param f                A generalized one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the one-electron operator.
     */
    template <typename Scalar>
    SquareMatrix<Scalar> evaluateOperatorDense(const ScalarGSQOneElectronOperator<Scalar>& f) const {

        if (f.numberOfOrbitals() != this->numberOfOrbitals()) {
            throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorDense(const ScalarGSQOneElectronOperator<double>&): The number of orbitals of this ONV basis and the given one-electron operator are incompatible.");
        }

        // Initialize a container for the dense matrix representation, and fill it with the general evaluation function.
        MatrixRepresentationEvaluationContainer<SquareMatrix<Scalar>> container {this->dimension()};
        this->evaluate<SquareMatrix<Scalar>>(f, container);

        return container.evaluation();
    }


    /**
     *  Calculate the dense matrix representation of a generalized two-electron operator in this ONV basis.
     *
     *  @param g                A generalized two-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the two-electron operator.
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarGSQTwoElectronOperator<double>& g) const;

    /**
     *  Calculate the dense matrix representation of a generalized Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A generalized Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the Hamiltonian.
     */
    SquareMatrix<double> evaluateOperatorDense(const GSQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Dense unrestricted operator evaluations
     */

    /**
     *  Calculate the dense matrix representation of a component of an unrestricted one-electron operator in this ONV basis.
     *
     *  @param f                A component of an unrestricted one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the one-electron operator.
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarUSQOneElectronOperatorComponent<double>& f) const;

    /**
     *  Calculate the dense matrix representation of a component of an unrestricted two-electron operator in this ONV basis.
     *
     *  @param g                A component of an unrestricted two-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the one-electron operator.
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarPureUSQTwoElectronOperatorComponent<double>& g) const;


    /*
     *  MARK: Diagonal generalized operator evaluations
     */

    /**
     *  Calculate the diagonal of the matrix representation of a generalized one-electron operator in this ONV basis.
     *
     *  @param f_op             A generalized one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the one-electron operator.
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarGSQOneElectronOperator<double>& f_op) const;

    /**
     *  Calculate the diagonal of the matrix representation of a generalized two-electron operator in this ONV basis.
     *
     *  @param g_op             A generalized two-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the two-electron operator.
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarGSQTwoElectronOperator<double>& g_op) const;

    /**
     *  Calculate the diagonal of the dense matrix representation of a generalized Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A generalized Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the Hamiltonian.
     */
    VectorX<double> evaluateOperatorDiagonal(const GSQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Sparse generalized operator evaluations
     */

    /**
     *  Calculate the sparse matrix representation of a generalized one-electron operator in this ONV basis.
     *
     *  @param f                A generalized one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A sparse matrix represention of the one-electron operator.
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarGSQOneElectronOperator<double>& f) const;

    /**
     *  Calculate the sparse matrix representation of a generalized two-electron operator in this ONV basis.
     *
     *  @param g                A generalized two-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A sparse matrix represention of the two-electron operator.
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarGSQTwoElectronOperator<double>& g) const;

    /**
     *  Calculate the sparse matrix representation of a generalized Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A generalized Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return A sparse matrix represention of the Hamiltonian.
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const GSQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Sparse unrestricted operator evaluations
     */

    /**
     *  Calculate the sparse matrix representation of a component of an unrestricted one-electron operator in this ONV basis.
     *
     *  @param f                A component of an unrestricted one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A sparse matrix represention of the one-electron operator.
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarUSQOneElectronOperatorComponent<double>& f) const;

    /**
     *  Calculate the sparse matrix representation of a component of an unrestricted two-electron operator in this ONV basis.
     *
     *  @param g                A component of an unrestricted two-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A sparse matrix represention of the one-electron operator.
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarPureUSQTwoElectronOperatorComponent<double>& g) const;


    /*
     *  MARK: Generalized matrix-vector product evaluations
     */

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a generalized one-electron operator with the given coefficient vector.
     *
     *  @param f                A generalized one-electron operator expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the one-electron operator.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarGSQOneElectronOperator<double>& f, const VectorX<double>& x) const;

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a generalized two-electron operator with the given coefficient vector.
     *
     *  @param g                A generalized two-electron operator expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the two-electron operator.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarGSQTwoElectronOperator<double>& g, const VectorX<double>& x) const;

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a generalized Hamiltonian with the given coefficient vector.
     *
     *  @param hamiltonian      A generalized Hamiltonian expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const GSQHamiltonian<double>& hamiltonian, const VectorX<double>& x) const;


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
        const auto N = this->numberOfElectrons();

        // Iterate over all ONVs, start with ONV with address 0.
        SpinUnresolvedONV onv = this->constructONVFromAddress(0);
        for (; !container.isFinished(); container.increment()) {
            for (size_t e = 0; e < N; e++) {  // Loop over all electrons that can be annihilated.

                auto q = onv.occupationIndexOf(e);  // Retrieve the orbital index of the electron that will be annihilated.

                // The diagonal values are a result of annihilation-creation on the same orbital index and are thus the same as the initial ONV.
                container.addRowwise(container.index, f(q, q));  // F(I,I)


                // For the non-diagonal values, we will create all possible matrix elements of the Hamiltonian in the routine below. Initialize an ONVPath that corresponds to the current ONV.
                ONVPath<SpinUnresolvedONVBasis> onv_path {*this, onv};
                onv_path.annihilate(q, e);

                // Stop the loop if:
                //      - 1) the path is finished, meaning that orbital index p is at M (the total number of orbitals)
                //      - and 2) if the orbital index is out of bounds after left translation of a vertical arc.
                while (!onv_path.isFinished() && onv_path.isOrbitalIndexValid()) {

                    // Find the next unoccupied orbital, i.e. the next vertical arc in the path.
                    onv_path.leftTranslateDiagonalArcUntilVerticalArc();

                    // Calculate the address of the path if we would close it right now.
                    const auto address = onv_path.addressAfterCreation();
                    const auto value = static_cast<double>(onv_path.sign()) * f(onv_path.orbitalIndex(), q);
                    const auto conjugated_value = GQCP::conj(value);  // For real numbers, this does nothing.

                    // Add the one-electron integral as matrix elements of a Hermitian matrix.
                    container.addRowwise(address, value);                // F(I,J)
                    container.addColumnwise(address, conjugated_value);  // F(J,I)

                    // Move the orbital index such that other unoccupied orbitals can be found within the loop.
                    onv_path.leftTranslateVerticalArc();
                }
            }

            // Prevent the last ONV since there is no possibility for an electron to be annihilated anymore.
            if (container.index < dim - 1) {
                this->transformONVToNextPermutation(onv);
            }
        }
    }


    /**
     *  Calculate the matrix representation of a generalized Hamiltonian in this ONV basis and emplace it in the given container.
     *
     *  @tparam Matrix                      The type of matrix used to store the evaluations.
     *
     *  @param hamiltonian                  An generalized Hamiltonian expressed in an orthonormal spinor basis.
     *  @param container                    A specialized container for emplacing evaluations/matrix elements.
     */
    template <typename Matrix>
    void evaluate(const GSQHamiltonian<double>& hamiltonian, MatrixRepresentationEvaluationContainer<Matrix>& container) const {

        // Prepare some variables.
        const auto& h_op = hamiltonian.core();
        const auto& g_op = hamiltonian.twoElectron();  // 'op' for 'operator'

        // For the most efficient evaluation of the one- and two-electron contributions combined, we use the effective one-electron operator.
        const auto k_op = g_op.effectiveOneElectronPartition() + h_op;
        const auto& k = k_op.parameters();

        const auto& g = g_op.parameters();


        const size_t M = this->numberOfOrbitals();
        const size_t N = this->numberOfElectrons();
        const size_t dim = this->dimension();


        SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0
        for (; !container.isFinished(); container.increment()) {   // I loops over all addresses in the spin-unresolved ONV basis
            if (container.index > 0) {
                this->transformONVToNextPermutation(onv);
            }
            int sign1 = -1;                      // start with -1 because we flip at the start of the annihilation (so we start at 1, followed by:  -1, 1, ...)
            for (size_t e1 = 0; e1 < N; e1++) {  // A1 (annihilation 1)

                sign1 *= -1;
                size_t p = onv.occupationIndexOf(e1);  // retrieve the index of a given electron
                size_t address = container.index - this->vertexWeight(p, e1 + 1);

                // Strictly diagonal values
                container.addRowwise(container.index, k(p, p));
                for (size_t q = 0; q < M; q++) {  // q loops over SOs
                    if (onv.isOccupied(q)) {
                        container.addRowwise(container.index, 0.5 * g(p, p, q, q));
                    } else {
                        container.addRowwise(container.index, 0.5 * g(p, q, q, p));
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

                        while (s < M) {
                            size_t J = address3 + this->vertexWeight(s, e4);
                            int signev = sign1 * sign2 * sign3 * sign4;
                            double value = signev * 0.5 * (g(p, q, r, s) + g(r, s, p, q) - g(p, s, r, q) - g(r, q, p, s));


                            container.addColumnwise(J, value);
                            container.addRowwise(J, value);

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
                while (q < M) {

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

                        while (s < M) {
                            size_t J = address3 + this->vertexWeight(s, e4);
                            int signev = sign1 * sign2 * sign3 * sign4;

                            double value = signev * 0.5 * (g(p, q, r, s) + g(r, s, p, q) - g(r, q, p, s) - g(p, s, r, q));

                            container.addColumnwise(J, value);
                            container.addRowwise(J, value);

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
                        while (s < M) {

                            size_t J = address2 + this->vertexWeight(s, e4);

                            int signev = sign1 * sign2 * sign3 * sign4;
                            double value = signev * 0.5 * (g(p, q, r, s) + g(r, s, p, q) - g(r, q, p, s) - g(p, s, r, q));

                            container.addColumnwise(J, value);
                            container.addRowwise(J, value);

                            s++;
                            this->shiftUntilNextUnoccupiedOrbital<1>(onv, address2, s, e4, sign4);
                        }
                    }

                    /**
                     *  A2 = C1
                     */
                    int signev = sign2 * sign1;

                    double value_I = k(p, q);  // cover the one electron calculations

                    for (size_t s = 0; s < M; s++) {
                        if (!onv.isOccupied(s)) {
                            value_I += 0.5 * (g(p, s, s, q));
                        } else {
                            value_I += 0.5 * (g(s, s, p, q) - g(s, q, p, s) + g(p, q, s, s));
                        }
                    }

                    value_I *= signev;

                    q++;

                    container.addColumnwise(address1, value_I);
                    container.addRowwise(address1, value_I);

                    this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign2);
                }
            }
        }
    }


    /*
     *  MARK: Legacy code
     */

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
};


}  // namespace GQCP
