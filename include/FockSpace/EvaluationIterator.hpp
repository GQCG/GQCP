// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#pragma once


#include "Mathematical/Representation/SquareMatrix.hpp"

#include <Eigen/Sparse>


namespace GQCP {


/**
 *  A templated private class for Fock spaces (private constructors with Fock spaces friends)
 *  It supports efficient dense or sparse storage for elements evaluated in the Fock space.
 *
 *  @tparam Matrix              the type of matrix in which the evaluations of the Fock space will be stored
 */
template<class Matrix>
class EvaluationIterator {
    size_t index = 0;
    size_t end;
    Matrix matrix;  // matrix containing the evaluations

    // CONSTRUCTOR
    /**
     * @param dimension         the dimensions of the matrix (equal to that of the fock space)
     */
    EvaluationMatrix(size_t dimension) :  matrix(Matrix::Zero(dimension, dimension)), end(dimension) {}


    // PUBLIC METHODS
    /**
     *  General interface to perform additions of values to elements of the matrix
     *
     * @param i         row index of the matrix
     * @param j         column index of the matrix
     * @param value     the value which is added to a given position in the matrix
     */
    void add_columnwise(size_t j, double value) {
        this->matrix(index, j) += value;
    }

    /**
     *  Method to signify the fact that the value is added to the lower quadrant of the given matrix evaluation.
     *  It exists to a satisfy the API used in the evaluation methods of Fock spaces.  This method exists to support matvec operations.
     *
     *  @param i         row index of the matrix
     *  @param j         column index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void add_rowwise(size_t j, double value) {
        this->matrix(j, index) += value;
    }

    // GETTER
    const Matrix& evaluation() const { return this->matrix; }

    // Friend Classes
    friend class FockSpace;
    friend class SelectedFockSpace;
    friend class ProductFockSpace;
};



/**
 *  Sparse template specialization is required because insertions into an existing sparse matrix are expensive
 *  Elements should only be added to the matrix once all of them are evaluated in a vector of triplets
 *  Therefore the "add(size_t, size_t, double)" method adds elements to a vector of triplets instead.
 */
template<>
class EvaluationMatrix<Eigen::SparseMatrix<double>> {

    Eigen::SparseMatrix<double> matrix;  // matrix containing the evaluations
    std::vector<Eigen::Triplet<double>> triplet_vector;  // vector which temporarily contains the added values

    // CONSTRUCTOR
    /**
     * @param dimension         the dimensions of the matrix (equal to that of the fock space)
     */
    EvaluationMatrix(size_t dimension) : matrix(Eigen::SparseMatrix<double>(dimension, dimension)) {}


    // PUBLIC METHODS
    /**
     *  Reserves an amount of memory for the triplet vector
     *
     *  @param n        amount of memory reserved
     */
    void reserve(size_t n) {
        this->triplet_vector.reserve(n);
        this->matrix.reserve(n);
    }

    /**
     *  General interface to perform additions of values to elements of the matrix
     *  This function adds the values to a triplet vector
     *  to add the values to the matrix, one should call "addToMatrix()"
     *
     *  @param i         row index of the matrix
     *  @param j         column index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void add(size_t i, size_t j, double value) {
        this->triplet_vector.emplace_back(i, j, value);
    }


    /**
     *  Method to signify the fact that the value is added to the lower quadrant of the given matrix evaluation
     *  It exists to a satisfy the API used in the evaluation methods of Fock spaces.  This method exists to support matvec operations.
     *
     *  @param i         row index of the matrix
     *  @param j         column index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void add_lower(size_t i, size_t j, double value) {
        this->triplet_vector.emplace_back(i, j, value);
    }


    /**
     *  Fill the sparse matrix with the elements stored in the triplet vector
     *  the more elements that are already present in the matrix, the more expensive this operation becomes
     *  Therefore, it is ill-advised to call the method more than once
     *
     *  After filling, the triplet vector is cleared
     */
    void addToMatrix() {
        this->matrix.setFromTriplets(this->triplet_vector.begin(), this->triplet_vector.end());
        this->triplet_vector = {};
    }

    /*
     *  This method does nothing, it exists to a satisfy the API used in the evaluation methods of Fock spaces.  This method exists to support matvec operations
     */   
    void flush(size_t index) {}


    // GETTERS
    const Eigen::SparseMatrix<double>& get_matrix() const { return this->matrix; }
    const std::vector<Eigen::Triplet<double>>& get_triplets() const { return triplet_vector; }


    // Friend Classes
    friend class FockSpace;
    friend class SelectedFockSpace;
    friend class ProductFockSpace;
};



/**
 *  Vector template specialization is required because matvec additions into an existing sparse matrix are expensive
 *  Elements should only be added to the matrix once all of them are evaluated in a vector of triplets
 *  Therefore the "add(size_t, size_t, double)" method adds elements to a vector of triplets instead.
 */
template<>
class EvaluationMatrix<VectorX<double>> {

    VectorX<double> matvec;  // matvec containing the evaluations
    const VectorX<double>& coefficient_vector;  // matvec containing the evaluations
    double sequential_double = 0;  // double which temporarily contains the added values, which can be flushed
    double nonsequential_double = 0;  // double gathered from the coefficient for nonsequential matvec additions.
    size_t index = 0;

    // CONSTRUCTOR
    /**
     * @param dimension         the dimensions of the matrix (equal to that of the fock space)
     */
    EvaluationMatrix(const VectorX<double>& coefficient_vector, const VectorX<double>& diagonal) : 
    coefficient_vector(coefficient_vector), 
    matvec(diagonal.cwiseProduct(coefficient_vector))
    {}

    EvaluationMatrix(const VectorX<double>& coefficient_vector) : 
    coefficient_vector(coefficient_vector), 
    matvec(VectorX<double>::Zero(coefficient_vector.rows()))
    {}

    // PUBLIC METHODS

    /**
     *  General interface to perform additions of values to elements of the matrix
     *  This function adds the values to a triplet vector
     *  to add the values to the matrix, one should call "addToMatrix()"
     *
     *  @param i         row index of the matrix
     *  @param j         column index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void add(size_t i, size_t j, double value) {
        sequential_double += value * coefficient_vector(j);
    }

    /**
     *  Method to signify the fact that the value is added to the lower quadrant of the given matrix evaluation
     *  It exists to a satisfy the API used in the evaluation methods of Fock spaces.  This method exists to support matvec operations.
     *
     *  @param i         row index of the matrix
     *  @param j         column index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void add_lower(size_t i, size_t j, double value) {
        this->matvec(i) += value * nonsequential_double;
    }

    /*
     *  Prepare the matvec variables for the iteration over the given index
     */  
    void flush(size_t index) {
        this->matvec(this->index) += sequential_double;
        this->index = index;
        sequential_double = 0;
        nonsequential_double = coefficient_vector(index);
    }


    // GETTERS
    const VectorX<double>& get_matrix() const { return this->matvec; }


    // Friend Classes
    friend class FockSpace;
    friend class SelectedFockSpace;
    friend class ProductFockSpace;
};

}  // namespace GQCP
