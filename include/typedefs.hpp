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
#ifndef GQCP_COMMON_HPP
#define GQCP_COMMON_HPP


#include <complex>
#include <cstdlib>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


/**
 *  A header that contains all typedefs inside the GQCP namespace
 */


namespace GQCP {


/*
 *  SCALARS
 */
using cd = std::complex<double>;


/*
 *  VECTORS
 */
using Vectoru = std::vector<size_t>;

template<typename Scalar>
using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

using VectorXs = Vector<size_t>;


/*
 *  MATRICES
 */
using Matrixu = std::vector<Vectoru>;

template<typename Scalar>
using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;


/*
 *  FUNCTIONS
 */
using VectorFunction = std::function<Eigen::VectorXd (const Eigen::VectorXd&)>;
using MatrixFunction = std::function<Eigen::MatrixXd (const Eigen::VectorXd&)>;


/*
 *  TENSORS
 */
template<typename Scalar>
using FourIndexTensor = Eigen::Tensor<Scalar, 4>;


/*
 *  FEATURES
 */
template<bool B, class T = void>
using enable_if_t = typename std::enable_if<B,T>::type;  // only implemented in C++14


}  // namespace GQCP


#endif  // GQCP_COMMON_HPP
