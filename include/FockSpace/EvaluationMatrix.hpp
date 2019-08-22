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
class EvaluationMatrix {

    Matrix matrix;  // matrix containing the evaluations

    // CONSTRUCTOR
    /**
     * @param dimension         the dimensions of the matrix (equal to that of the fock space)
     */
    EvaluationMatrix(size_t dimension) :  matrix(Matrix::Zero(dimension, dimension)) {}


    // PUBLIC METHODS
    /**
     *  General interface to perform additions of values to elements of the matrix
     *
     * @param i         row index of the matrix
     * @param j         column index of the matrix
     * @param value     the value which is added to a given position in the matrix
     */
    void add(size_t i, size_t j, double value) {
        this->matrix(i, j) += value;
    }

    // GETTER
    const Matrix& get_matrix() const { return this->matrix; }

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


    // GETTERS
    const Eigen::SparseMatrix<double>& get_matrix() const { return this->matrix; }
    const std::vector<Eigen::Triplet<double>>& get_triplets() const { return triplet_vector; }


    // Friend Classes
    friend class FockSpace;
    friend class SelectedFockSpace;
    friend class ProductFockSpace;
};


}  // namespace GQCP
