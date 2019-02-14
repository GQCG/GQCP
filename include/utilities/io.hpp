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
#ifndef GQCP_IO_HPP
#define GQCP_IO_HPP


#include <string>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>



namespace GQCP {


/**
 *  Read a vector from a given file and put the elements in the given vector
 *
 *  @param filename     the name of the file to be read in
 *  @param v            the vector that gets the contents of the file
 */
void readVectorFromFile(const std::string& filename, Eigen::VectorXd& v);

/**
 *  Read a matrix from a given file and put the elements in the given matrix
 *
 *  @param filename     the name of the file to be read in
 *  @param M            the matrix that gets the contents of the file
 */
void readArrayFromFile(const std::string& filename, Eigen::MatrixXd& M);

/**
 *  Read a tensor from a given file and put the elements in the given tensor
 *
 *  @param filename     the name of the file to be read in
 *  @param T            the tensor that gets the contents of the file
 */
void readArrayFromFile(const std::string& filename, Eigen::Tensor<double, 4>& T);

/**
 *  @param M        the matrix to be printed
 *
 *  Print the contents of a matrix in a fashionable way
 */
void print(const Eigen::MatrixXd& M);

/**
 *  Print the contents of a matrix to an output filestream
 *
 *  @param M                        the matrix to be printed
 *  @param output_filestream        the stream used for outputting
 */
void print(const Eigen::MatrixXd& M, std::ofstream& output_filestream);

/**
 *  Print the contents of a matrix to an output stream
 *
 *  @param T                        the tensor to be printed
 *  @param output_stream            the stream used for outputting
 */
void print(const Eigen::Tensor<double, 4>& T, std::ostream& output_stream);

/**
 *  Print the contents of a tensor in a fashionable way
 *
 *  @param T        the tensor to be printed
 */
void print(const Eigen::Tensor<double, 4>& T);


}  // namespace GQCP



#endif  // GQCP_IO_HPP
