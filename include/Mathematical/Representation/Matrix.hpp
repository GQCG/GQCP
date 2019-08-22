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


#include "Mathematical/CartesianDirection.hpp"
#include "typedefs.hpp"

#include <Eigen/Dense>

#include <boost/algorithm/string.hpp>

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
class Matrix : public Eigen::Matrix<_Scalar, _Rows, _Cols> {
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
     *  Construct a vector by reading in a file
     *
     *  @param filename     the name of the file to be read in
     *  @param rows         the number of expected rows
     */
    template <typename Z = Self>  // enable_if must have Z inside
    static enable_if_t<Self::is_vector, Z> FromFile(const std::string& filename, size_t rows) {

        // Initialize a zero vector
        Self result = Self::Zero(rows);

        std::ifstream file (filename);
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

        std::ifstream file (filename);
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

                result(i,j) = value;
            }

            file.close();
        } else {
            throw std::runtime_error("Matrix::FromFile(std::string, size_t, size_t): Cannot open the given file. Maybe you specified a wrong path?");
        }

        return result;
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
    enable_if_t<Self::is_vector && (Rows == 3), Z> operator()(CartesianDirection direction) const {
        
        const auto& index = static_cast<size_t>(direction);  // 0, 1, or 2
        return this->operator()(index);
    }

    /**
     *  @param direction            the Cartesian direction (x, y, or z)
     * 
     *  @return a modifiable value in the vector that corresponds to the given direction
     */
    template <typename Z = Scalar&>
    enable_if_t<Self::is_vector && (Rows == 3), Z> operator()(CartesianDirection direction) {
        
        const auto& index = static_cast<size_t>(direction);  // 0, 1, or 2
        return this->operator()(index);
    }

    using Base::operator();  // bring over the other operator() overloads



    /*
     *  PUBLIC METHODS
     */


    /**
     *  @return this as a const Eigen::Matrix
     */
    const Base& Eigen() const {
        return static_cast<const Base&>(*this);
    }


    /**
     *  @return this as a non-const Eigen::Matrix
     */
    Base& Eigen() {
        return static_cast<Base&>(*this);
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
                output_stream << i << ' ' << j << "  " << this->operator()(i,j) << std::endl;
            }
        }
    }


    /**
     *  @param i    row index (starting from 0)
     *  @param j    column index (starting from 0)
     *
     *  @return the i-j minor (i.e. delete the i-th row and j-th column)
     */
    template <typename Z = Self>  // enable_if must have Z inside
    enable_if_t<Self::is_matrix, Z> matrixMinor(size_t i, size_t j) const {  // wrap minor inside braces as a fix for the GCC macro 'minor'

        // Delete the i-th row
        Self A_i = Self::Zero(this->rows() - 1, this->cols());

        for (size_t i2 = 0; i2 < this->rows(); i2++) {
            if (i2 < i) {
                A_i.row(i2) = this->row(i2);
            } else if (i2 == i) {
                continue;
            } else if (i2 > i) {
                A_i.row(i2-1) = this->row(i2);
            }
        }

        // Delete the j-th column
        Self A_ij = Self::Zero(this->rows() - 1, this->cols() - 1);
        for (size_t j2 = 0; j2 < this->cols(); j2++) {
            if (j2 < j) {
                A_ij.col(j2) = A_i.col(j2);
            } else if (j2 == j) {
                continue;
            } else if (j2 > j) {
                A_ij.col(j2-1) = A_i.col(j2);
            }
        }

        return A_ij;
    }



};



/*
 *  Convenience typedefs related to Matrix
 */

template <typename Scalar>
using MatrixX = Matrix<Scalar, Dynamic, Dynamic>;

template<typename Scalar, int Rows>
using Vector = Matrix<Scalar, Rows, 1>;

template <typename Scalar>
using VectorX = Vector<Scalar, Dynamic>;

using VectorXs = VectorX<size_t>;


using VectorFunction = std::function<VectorX<double> (const VectorX<double>&)>;
using MatrixFunction = std::function<MatrixX<double> (const VectorX<double>&)>;


}  // namespace GQCP
