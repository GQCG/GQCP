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


#include "Mathematical/CartesianDirection.hpp"
#include "Utilities/type_traits.hpp"
#include "Utilities/typedefs.hpp"

#include <boost/algorithm/string.hpp>

#include <Eigen/Dense>

#include <fstream>
#include <iostream>


namespace GQCP {


constexpr auto Dynamic = Eigen::Dynamic;

/**
 *  An extension of the Eigen::Matrix class, with extra operations
 *
 *  @tparam _Scalar     the scalar representation type
 *  @tparam _Rows       the number of rows (int or Dynamic)
 *  @tparam _Cols       the number of columns (int or Dynamic)
 *
 *  We have decided to inherit from Eigen::Matrix, because we will use different hierarchies: see also: https://eigen.tuxfamily.org/dox-devel/TopicCustomizing_InheritingMatrix.html
 */
template <typename _Scalar = double, int _Rows = Dynamic, int _Cols = Dynamic>
class Matrix: public Eigen::Matrix<_Scalar, _Rows, _Cols> {
public:
    using Scalar = _Scalar;
    static constexpr auto Rows = _Rows;
    static constexpr auto Cols = _Cols;

    using Self = Matrix<Scalar, Rows, Cols>;
    using Base = Eigen::Matrix<Scalar, Rows, Cols>;


private:
    static constexpr bool is_vector = (Cols == 1);
    static constexpr bool is_matrix = (Cols >= 2) || (Cols == Dynamic);


public:
    /*
     *  CONSTRUCTORS
     */

    using Eigen::Matrix<Scalar, _Rows, _Cols>::Matrix;  // inherit base constructors


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Convert a given column-major vector to a matrix with the given number of rows
     * 
     *  @param v            a vector that is supposed to be in a column-major ordering
     *  @param rows         the number of rows the resulting matrix should have
     *  @param cols         the number of columns the resulting matrix should have
     */
    template <typename Z = Self>
    static enable_if_t<(Cols == Dynamic) && (Rows == Dynamic), Z> FromColumnMajorVector(const Matrix<Scalar, Dynamic, 1>& v, const size_t rows, const size_t cols) {

        return Self(Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(v.data(), rows, cols));
    }


    /**
     *  Construct a vector by reading in a file
     *
     *  @param filename     the name of the file to be read in
     *  @param rows         the number of expected rows
     */
    template <typename Z = Self>  // enable_if must have Z inside
    static enable_if_t<Self::is_vector, Z> FromFile(const std::string& filename, size_t rows) {

        // Initialize a zero vector
        Self result = Self::Zero(rows);

        std::ifstream file {filename};
        size_t index = 0;
        if (file.is_open()) {
            std::string line;

            while (std::getline(file, line)) {
                std::vector<std::string> splitted_line;  // create a container for the line to be split in

                // Split the line on any whitespace or tabs.
                boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

                if (splitted_line.size() != 1) {
                    throw std::runtime_error("Matrix::FromFile(std::string, size_t): Found a line that doesn't contain exactly 1 field delimited by whitespace.");
                }

                auto value = std::stod(splitted_line[0]);
                result(index) = value;

                ++index;
            }

            file.close();
        } else {
            throw std::runtime_error("Matrix::FromFile(std::string, size_t): Cannot open the given file. Maybe you specified a wrong path?");
        }

        return result;
    }


    /**
     *  Construct a matrix by reading in a file
     *
     *  @param filename     the name of the file to be read in
     *  @param rows         the number of expected rows
     *  @param cols         the number of expected columns
     */
    template <typename Z = Self>  // enable_if must have Z inside
    static enable_if_t<Self::is_matrix, Z> FromFile(const std::string& filename, size_t rows, size_t cols) {

        Self result = Self::Zero(rows, cols);

        std::ifstream file {filename};
        if (file.is_open()) {
            std::string line;

            while (std::getline(file, line)) {
                std::vector<std::string> splitted_line;  // create a container for the line to be split in

                // Split the line on any whitespace or tabs.
                boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

                if (splitted_line.size() != 3) {
                    throw std::runtime_error("Matrix::FromFile(std::string, size_t, size_t): Found a line that doesn't contain exactly 3 fields delimited by whitespace.");
                }

                auto i = std::stoi(splitted_line[0]);
                auto j = std::stoi(splitted_line[1]);
                auto value = std::stod(splitted_line[2]);

                result(i, j) = value;
            }

            file.close();
        } else {
            throw std::runtime_error("Matrix::FromFile(std::string, size_t, size_t): Cannot open the given file. Maybe you specified a wrong path?");
        }

        return result;
    }


    /**
     *  Convert a given row-major vector to a matrix with the given number of rows
     * 
     *  @param v            a vector that is supposed to be in a column-major ordering
     *  @param rows         the number of rows the resulting matrix should have
     *  @param cols         the number of columns the resulting matrix should have
     */
    template <typename Z = Self>
    static enable_if_t<(Cols == Dynamic) && (Rows == Dynamic), Z> FromRowMajorVector(const Matrix<Scalar, Dynamic, 1>& v, const size_t rows, const size_t cols) {

        // Change the given vector's elements to a column-major representation
        Matrix<Scalar, Dynamic, 1> v_column_major = Matrix<Scalar, Dynamic, 1>::Zero(rows * cols);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                size_t row_major_index = i * cols + j;
                size_t column_major_index = j * rows + i;

                v_column_major(column_major_index) = v(row_major_index);
            }
        }

        return Self::FromColumnMajorVector(v_column_major, rows, cols);
    }


    /*
     *  OPERATORS
     */

    /**
     *  @param direction            the Cartesian direction (x, y, or z)
     * 
     *  @return the value in the vector that corresponds to the given direction
     */
    template <typename Z = Scalar>
    enable_if_t<Self::is_vector && (Rows == 3), Z> operator()(const CartesianDirection direction) const {

        const auto index = static_cast<size_t>(direction);  // 0, 1, or 2
        return this->operator()(index);
    }

    /**
     *  @param direction            the Cartesian direction (x, y, or z)
     * 
     *  @return a modifiable value in the vector that corresponds to the given direction
     */
    template <typename Z = Scalar&>
    enable_if_t<Self::is_vector && (Rows == 3), Z> operator()(const CartesianDirection direction) {

        const auto index = static_cast<size_t>(direction);  // 0, 1, or 2
        return this->operator()(index);
    }

    using Base::operator();  // bring over the other operator() overloads


    /*
     *  PUBLIC METHODS
     */


    /**
     *  @param i    row index (starting from 0)
     *  @param j    column index (starting from 0)
     *
     *  @return the i-j minor (i.e. delete the i-th row and j-th column)
     */
    template <typename Z = Self>
    enable_if_t<Self::is_matrix, Z> calculateMinor(size_t i, size_t j) const {

        Self this_copy = *this;

        this_copy.removeRow(i);
        this_copy.removeColumn(j);

        return this_copy;
    }


    /**
     *  @return this as a const Eigen::Matrix
     */
    const Base& Eigen() const { return static_cast<const Base&>(*this); }

    /**
     *  @return this as a non-const Eigen::Matrix
     */
    Base& Eigen() { return static_cast<Base&>(*this); }


    /**
     *  @param start_i      the index at which the rows should start
     *  @param start_j      the index at which the columns should start
     * 
     *  @return this matrix as a vector in a column-major storage, i.e. a pair-wise reduced form of this matrix. The elements of the matrix are put into the vector such that
     *      v(m) = M(i,j)
     *
     *  in which
     *      m is calculated from i and j in a column-major way
     */
    template <typename Z = Matrix<Scalar, Rows, 1>>
    enable_if_t<Self::is_matrix, Z> pairWiseReduced(const size_t start_i = 0, size_t start_j = 0) const {

        // Initialize the resulting vector
        Matrix<Scalar, Rows, 1> v((this->rows() - start_i) * (this->cols() - start_j));

        // Calculate the compound indices and bring the elements from the matrix over into the vector
        size_t vector_index = 0;
        for (size_t j = start_j; j < this->cols(); j++) {      // "column major" ordering for row_index<-i,j so we do j first, then i
            for (size_t i = start_i; i < this->rows(); i++) {  // in column major indices, columns are contiguous, so the first of two indices changes more rapidly

                v(vector_index) = this->operator()(i, j);
                vector_index++;
            }
        }

        return v;
    }


    /**
     *  Print the contents of a this to an output filestream
     *
     *  @param output_stream        the stream used for outputting
     */
    template <typename Z = void>  // enable_if must have Z inside
    enable_if_t<Self::is_matrix, Z> print(std::ostream& output_stream = std::cout) const {

        for (size_t i = 0; i < this->rows(); i++) {
            for (size_t j = 0; j < this->cols(); j++) {
                output_stream << i << ' ' << j << "  " << this->operator()(i, j) << std::endl;
            }
        }
    }


    /**
     *  Remove the i-th column of this matrix.
     * 
     *  @param i            the index of the column that should be removed
     * 
     *  @note Implementation adapted from (https://stackoverflow.com/a/46303314), while waiting for Eigen's 3.4 release.
     */
    void removeColumn(const size_t i) {

        const size_t rows = this->rows();         // the number of rows of the resulting matrix
        const size_t columns = this->cols() - 1;  // the number of columns of the resulting matrix

        // Shift the columns to the left.
        if (i < columns) {
            this->block(0, i, rows, columns - i) = this->rightCols(columns - i);
        }

        this->conservativeResize(rows, columns);
    }


    /**
     *  Remove the columns at the given indices.
     * 
     *  @param column_indices               the indices of the columns that should be removed, in ascending order
     */
    void removeColumns(const std::vector<size_t>& column_indices) {

        size_t removed_counter = 0;
        for (const auto& i : column_indices) {
            this->removeColumn(i - removed_counter);  // take into account that the given indices are with respect to the old dimensions of the matrix
            removed_counter++;
        }
    }


    /**
     *  Remove the i-th row of this matrix.
     * 
     *  @param i            the index of the row that should be removed
     * 
     *  @note Implementation adapted from (https://stackoverflow.com/a/46303314), while waiting for Eigen's 3.4 release.
     */
    void removeRow(const size_t i) {

        const size_t rows = this->rows() - 1;  // the number of rows of the resulting matrix
        const size_t columns = this->cols();   // the number of columns of the resulting matrix

        // Shift the rows up.
        if (i < rows) {
            this->block(i, 0, rows - i, columns) = this->bottomRows(rows - i);
        }

        this->conservativeResize(rows, columns);
    }


    /**
     *  Remove the rows at the given indices.
     * 
     *  @param row_indices                  the indices of the columns that should be removed, in ascending order
     */
    void removeRows(const std::vector<size_t>& row_indices) {

        size_t removed_counter = 0;
        for (const auto& i : row_indices) {
            this->removeRow(i - removed_counter);  // take into account that the given indices are with respect to the old dimensions of the matrix
            removed_counter++;
        }
    }
};


/*
 *  Convenience typedefs related to Matrix. These have to be in order.
 */

template <typename Scalar, int Rows>
using Vector = Matrix<Scalar, Rows, 1>;

template <typename Scalar>
using MatrixX = Matrix<Scalar, Dynamic, Dynamic>;

template <typename Scalar>
using VectorX = Vector<Scalar, Dynamic>;

template <typename Scalar>
using VectorFunction = std::function<VectorX<Scalar>(const VectorX<Scalar>&)>;

template <typename Scalar>
using MatrixFunction = std::function<MatrixX<Scalar>(const VectorX<Scalar>&)>;


}  // namespace GQCP
