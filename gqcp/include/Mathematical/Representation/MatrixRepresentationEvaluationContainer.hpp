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


#include "Mathematical/Representation/SquareMatrix.hpp"

#include <Eigen/Sparse>


namespace GQCP {


/**
 *  A class that can efficiently construct the matrix representation of a second-quantized operator.
 * 
 *  @tparam _Matrix             The type that is used as a matrix representation, i.e. a container that stores the matrix representation.
 */
template <typename _Matrix>
class MatrixRepresentationEvaluationContainer {
public:
    // The type that is used as a matrix representation, i.e. a container that stores the matrix representation.
    using Matrix = _Matrix;

    // The scalar type of one of the matrix elements: real or complex.
    using Scalar = typename Matrix::Scalar;


public:
    // The current index of the iterator, which is an address in a particular ONV basis.
    size_t index = 0;

    // The expected last position of the iterator, i.e. the dimension of the particular ONV basis.
    size_t end;

    // The matrix representation.
    Matrix matrix;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Initialize an empty container, ready to start an evaluation.
     * 
     *  @param dimension            The dimension of the matrix representation, i.e. the number of elements in one row/column.
     */
    MatrixRepresentationEvaluationContainer(const size_t dimension) :
        matrix {Matrix::Zero(dimension)},
        end {dimension} {}


    /*
     *  MARK: Additions and evaluations
     */

    /**
     *  Add a value to the diagonal of the matrix representation. The current iterator index corresponds to the position on the diagonal.
     * 
     *  @param value            The value.
     */
    void addDiagonally(const Scalar value) {
        this->matrix(this->index, this->index) += value;
    }

    /**
     *  Add a value to the matrix representation, where the current iterator index corresponds to the row.
     * 
     *  @param column           The column index of the matrix representation at which the given value should be added.
     *  @param value            The value.
     */
    void addColumnwise(const size_t column, const Scalar value) {
        this->matrix(this->index, column) += value;
    }


    /**
     *  Add a value to the matrix representation, where the current iterator index corresponds to the column.
     * 
     *  @param row              The row index of the matrix representation at which the given value should be added.
     *  @param value            The value.
     */
    void addRowwise(const size_t row, const Scalar value) {
        this->matrix(row, this->index) += value;
    }


    /**
     *  @return The matrix representation.
     */
    const Matrix& evaluation() const { return this->matrix; }


    /*
     *  MARK: Iterator
     */

    /**
     *  Move to the next index in the iteration.
     */
    void increment() { this->index++; }

    /**
     *  Check if the iteration is finished. If so, the iterator is prepared to start a new evaluation.
     * 
     *  @return If the iteration is finished.
     */
    bool isFinished() {

        if (this->index == this->end) {
            this->index = 0;
            return true;
        } else {
            return false;
        }
    }
};


/**
 *  Sparse template specialization is required because insertions into an existing sparse matrix are expensive
 *  Elements should only be added to the matrix once all of them are evaluated in a vector of triplets
 * 
 *  @tparam _Scalar         The scalar type of one of the matrix elements: real or complex.
 */
template <typename _Scalar>
class MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<_Scalar>> {
public:
    // The scalar type of one of the matrix elements: real or complex.
    using Scalar = _Scalar;


public:
    size_t index = 0;  // current position of the iterator in the dimension of the ONV basis
    size_t end;        // total dimension

    Eigen::SparseMatrix<Scalar> matrix;                  // matrix containing the evaluations
    std::vector<Eigen::Triplet<Scalar>> triplet_vector;  // vector which temporarily contains the added values


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     * @param dimension         the dimensions of the matrix (equal to that of the fock space)
     */
    MatrixRepresentationEvaluationContainer(const size_t dimension) :
        matrix {Eigen::SparseMatrix<Scalar>(dimension, dimension)},
        end(dimension) {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Add a value to the matrix evaluation in which the current iterator index corresponds to the row and the given index corresponds to the column
     *  This function adds the values to a triplet vector
     *  to add the values to the sparse matrix, one should call "addToMatrix()"
     * 
     *  @param column    column index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void addColumnwise(const size_t column, const Scalar value) { this->triplet_vector.emplace_back(this->index, column, value); }

    /**
     *  Add a value to the matrix evaluation in which the current iterator index corresponds to the column and the given index corresponds to the row
     *  This function adds the values to a triplet vector
     *  to add the values to the sparse matrix, one should call "addToMatrix()"
     * 
     *  @param row       row index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void addRowwise(const size_t row, const Scalar value) { this->triplet_vector.emplace_back(row, this->index, value); }

    /**
     *  Fill the sparse matrix with the elements stored in the triplet vector
     *  the more elements that are already present in the matrix, the more expensive this operation becomes
     *  Therefore, it is ill-advised to call the method more than once,
     *  after filling, the triplet vector is cleared
     */
    void addToMatrix() {
        this->matrix.setFromTriplets(this->triplet_vector.begin(), this->triplet_vector.end());
        this->triplet_vector = {};
    }

    /**
     *  @return the evaluation that is stored
     * 
     *  @note the matrix will not be initialized if `addToMatrix()` has not been called
     */
    const Eigen::SparseMatrix<Scalar>& evaluation() const { return this->matrix; }

    /**
     *  Move to the next index in the iteration
     */
    void increment() { this->index++; }

    /**
     *  Tests if the iteration is finished, if true the index is reset to 0 
     * 
     *  @return true if the iteration is finished
     */
    bool isFinished() {
        if (this->index == this->end) {
            this->index = 0;
            return true;
        } else {
            return false;
        }
    }

    /**
     *  Reserves an amount of memory for the triplet vector
     *
     *  @param n        the amount triplets that should be reserved
     */
    void reserve(const size_t n) {
        this->triplet_vector.reserve(n);
        this->matrix.reserve(n);
    }

    /**
     *  @return the triplet vector
     */
    const std::vector<Eigen::Triplet<Scalar>>& triplets() const { return triplet_vector; }
};


/**
 *  Vector template specialization is required because of matvec evaluations are stored in a vector additions
 * 
 *  @tparam _Scalar         The scalar type of one of the matrix elements: real or complex.
 */
template <typename _Scalar>
class MatrixRepresentationEvaluationContainer<VectorX<_Scalar>> {
public:
    // The scalar type of one of the matrix elements: real or complex.
    using Scalar = _Scalar;

public:
    size_t index = 0;  // current position of the iterator in the dimension of the ONV basis
    size_t end;        // total dimension

    VectorX<Scalar> matvec;                     // matvec containing the evaluations
    const VectorX<Scalar>& coefficient_vector;  // vector with which is multiplied
    double sequential_double = 0;               // double which temporarily contains the sum of added values, it corresponds to the value in the matvec of the current index and allows a for a single write access each iteration instead of multiple
    double nonsequential_double = 0;            // double gathered from the coefficient for non-sequential matvec additions, this corresponds to the value of the coefficient vector of the current index, it allows for a single read operation each iteration


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     * @param dimension         the dimensions of the matrix (equal to that of the fock space)
     */
    MatrixRepresentationEvaluationContainer(const VectorX<Scalar>& coefficient_vector, const VectorX<Scalar>& diagonal) :
        end {static_cast<size_t>(coefficient_vector.rows())},
        coefficient_vector {coefficient_vector},
        matvec {diagonal.cwiseProduct(coefficient_vector)} {}

    MatrixRepresentationEvaluationContainer(const VectorX<Scalar>& coefficient_vector) :
        end {static_cast<size_t>(coefficient_vector.rows())},
        coefficient_vector {coefficient_vector},
        matvec {VectorX<Scalar>::Zero(coefficient_vector.rows())} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Add a value to the matrix evaluation in which the current iterator index corresponds to the row and the given index corresponds to the column
     * 
     *  @param column    column index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void addColumnwise(const size_t column, const Scalar value) { sequential_double += value * coefficient_vector(column); }

    /**
     *  Add a value to the matrix evaluation in which the current iterator index corresponds to the column and the given index corresponds to the row
     * 
     *  @param row       row index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void addRowwise(const size_t row, const Scalar value) { this->matvec(row) += value * this->nonsequential_double; }

    /**
     *  @return the evaluation that is stored
     */
    const VectorX<Scalar>& evaluation() const { return this->matvec; }

    /**
     *  Move to the next index in the iteration, this is accompanied by an addition to the matvec and reset of the sequential double
     */
    void increment() {
        this->matvec(this->index) += this->sequential_double;
        this->sequential_double = 0;
        this->index++;
    }

    /**
     *  Tests if the iteration is finished, if true the index is reset to 0 
     *  If false the nonsequential_double is updated to the value of the current iteration
     * 
     *  @return true if the iteration is finished
     */
    bool isFinished() {
        if (this->index == this->end) {
            this->index = 0;
            return true;
        } else {
            this->nonsequential_double = this->coefficient_vector(this->index);
            return false;
        }
    }
};

}  // namespace GQCP
