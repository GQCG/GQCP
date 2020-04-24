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


#include "Mathematical/Representation/Matrix.hpp"
#include "Utilities/type_traits.hpp"

#include <unsupported/Eigen/CXX11/Tensor>


namespace GQCP {


/**
 *  An extension of the Eigen::Tensor class, with extra operations
 *
 *  @tparam _Scalar     the scalar representation type
 *  @tparam _Rank       the rank of the tensor, i.e. the number of indices
 *
 *  We have decided to inherit from Eigen::Tensor, because we will use different hierarchies: see also: https://eigen.tuxfamily.org/dox-devel/TopicCustomizing_InheritingMatrix.html
 */
template <typename _Scalar, int _Rank>
class Tensor: public Eigen::Tensor<_Scalar, _Rank> {
public:
    using Scalar = _Scalar;
    static constexpr auto Rank = _Rank;
    using Self = Tensor<Scalar, Rank>;
    using Base = Eigen::Tensor<Scalar, Rank>;


public:
    /*
     *  CONSTRUCTORS
     */

    using Eigen::Tensor<Scalar, _Rank>::Tensor;  // inherit base constructors


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  @param T            a rank-4 tensor
     *
     *  @param i            1st starting index
     *  @param j            2nd starting index
     *  @param k            3rd starting index
     *  @param l            4th starting index
     *  @param desize       early cut-off of index iteration
     *
     *  @return a rank-4 tensor from an other rank-4 tensor, starting from given indices
     */
    template <int Z = Rank>
    static enable_if_t<Z == 4, Self> FromBlock(const Self& T, size_t i, size_t j, size_t k, size_t l, size_t desize = 0) {

        Tensor<double, Rank> T_result {static_cast<long>(T.dimension(0) - i - desize),
                                       static_cast<long>(T.dimension(1) - j - desize),
                                       static_cast<long>(T.dimension(2) - k - desize),
                                       static_cast<long>(T.dimension(3) - l - desize)};
        T_result.setZero();

        for (size_t p = 0; p < T_result.dimension(0); p++) {
            for (size_t q = 0; q < T_result.dimension(1); q++) {
                for (size_t r = 0; r < T_result.dimension(2); r++) {
                    for (size_t s = 0; s < T_result.dimension(3); s++) {
                        T_result(p, q, r, s) = T(i + p, j + q, k + r, l + s);
                    }
                }
            }
        }

        return T_result;
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return this as a const Eigen::Tensor, as a work-around to fix Eigen::Tensor expressions
     */
    const Base& Eigen() const {
        return static_cast<const Base&>(*this);
    }

    /**
     *  @return this as a non-const Eigen::Tensor, as a work-around to fix Eigen::Tensor expressions
     */
    Base& Eigen() {
        return static_cast<Base&>(*this);
    }


    /**
     *  @param other        the other tensor
     *
     *  @return if this tensor has the same dimensions as the other tensor
     */
    bool hasEqualDimensions(const Self& other) const {

        for (size_t i = 0; i < Rank; i++) {
            if (this->dimension(i) != other.dimension(i)) {
                return false;
            }
        }

        return true;
    }


    /**
     *  @param other        the other tensor
     *  @param tolerance    the tolerance for element-wise comparison
     *
     *  @return if this is approximately equal to the other
     */
    template <int Z = Rank>
    enable_if_t<Z == 4, bool> isApprox(const Self& other, double tolerance = 1.0e-12) const {

        if (!this->hasEqualDimensions(other)) {
            throw std::invalid_argument("RankFourTensor<Scalar>::isApprox(Self, double): the tensors have different dimensions");
        }

        // Check every pair of values
        for (size_t i = 0; i < this->dimension(0); i++) {
            for (size_t j = 0; j < this->dimension(1); j++) {
                for (size_t k = 0; k < this->dimension(2); k++) {
                    for (size_t l = 0; l < this->dimension(3); l++) {
                        if (std::abs(this->operator()(i, j, k, l) - other(i, j, k, l)) > tolerance) {
                            return false;
                        }
                    }
                }
            }
        }  // rank-4 tensor traversing

        return true;
    }


    /**
     *  Print the contents of this to an output stream
     *
     *  @param output_stream        the stream used for outputting
     */
    template <int Z = Rank>
    enable_if_t<Z == 4> print(std::ostream& output_stream = std::cout) const {

        for (size_t i = 0; i < this->dimension(0); i++) {
            for (size_t j = 0; j < this->dimension(1); j++) {
                for (size_t k = 0; k < this->dimension(2); k++) {
                    for (size_t l = 0; l < this->dimension(3); l++) {
                        output_stream << i << ' ' << j << ' ' << k << ' ' << l << "  " << this->operator()(i, j, k, l) << std::endl;
                    }
                }
            }
        }
    }


    /**
     *  @param start_i      the index at which the first rank should start
     *  @param start_j      the index at which the second rank should start
     *  @param start_k      the index at which the third rank should start
     *  @param start_l      the index at which the fourth rank should start
     * 
     *  @return a pair-wise reduced form of this rank-4 tensor. The elements of the tensor are put into the matrix such that
     *      M(m,n) = T(i,j,k,l)
     *
     *  in which
     *      m is calculated from i and j in a column-major way
     *      n is calculated from k and l in a column-major way
     */
    template <int Z = Rank>
    enable_if_t<Z == 4, Matrix<Scalar>> pairWiseReduce(const size_t start_i = 0, const size_t start_j = 0, const size_t start_k = 0, const size_t start_l = 0) const {

        // Initialize the resulting matrix
        const auto dims = this->dimensions();
        Matrix<Scalar> M {(dims[0] - start_i) * (dims[1] - start_j),
                          (dims[2] - start_k) * (dims[3] - start_l)};

        // Calculate the compound indices and bring the elements from the tensor over into the matrix
        size_t row_index = 0;
        for (size_t j = start_j; j < dims[1]; j++) {      // "column major" ordering for row_index<-i,j so we do j first, then i
            for (size_t i = start_i; i < dims[0]; i++) {  // in column major indices, columns are contiguous, so the first of two indices changes more rapidly

                size_t column_index = 0;
                for (size_t l = start_l; l < dims[3]; l++) {      // "column major" ordering for column_index<-k,l so we do l first, then k
                    for (size_t k = start_k; k < dims[2]; k++) {  // in column major indices, columns are contiguous, so the first of two indices changes more rapidly

                        M(row_index, column_index) = this->operator()(i, j, k, l);

                        column_index++;
                    }
                }

                row_index++;
            }
        }

        return M;
    }


    /**
     *  Add a rank-4 tensor into this, starting from given indices
     *
     *  @param T        a rank-4 tensor
     *  @param i                starting index for the 1st index axis of the tensor
     *  @param j                starting index for the 2nd index axis of the tensor
     *  @param k                starting index for the 3rd index axis of the tensor
     *  @param l                starting index for the 4th index axis of the tensor
     *
     *  @return a reference to updated this
     */
    template <int Z = Rank>
    enable_if_t<Z == 4, Self&> addBlock(const Self& T, size_t i, size_t j, size_t k, size_t l) {

        for (size_t p = 0; p < T.dimension(0); p++) {
            for (size_t q = 0; q < T.dimension(1); q++) {
                for (size_t r = 0; r < T.dimension(2); r++) {
                    for (size_t s = 0; s < T.dimension(3); s++) {
                        this->operator()(i + p, j + q, k + r, l + s) += T(p, q, r, s);
                    }
                }
            }
        }

        return (*this);
    }


    /**
     *  Add a matrix to a this tensor starting from given indices
     *
     *  @tparam r               indicates with which tensor index axis (0,1,2,3) the row index axis of the matrix should align 
     *  @tparam s               indicates with which tensor index axis (0,1,2,3) the column index axis of the matrix should align 
     *
     *  @param M                a matrix
     *  @param i                starting index for the 1st index axis of the tensor
     *  @param j                starting index for the 2nd index axis of the tensor
     *  @param k                starting index for the 3rd index axis of the tensor
     *  @param l                starting index for the 4th index axis of the tensor
     *
     *  @return a reference to updated this
     *
     *
     *  Example:
     *      Given a rank-4 tensor of dimensions (10,10,10,10), and a matrix M of dimensions (3,3)
     *       Input : <2,0> (M, 0, 2, 1, 3):
     *       <2,0> dictates that the row index axis of the matrix aligns with the 3rd index axis of the tensor (2nd starting from 0)
     *       and that the column index axis of the matrix aligns with the 1st index axis tensor (0th starting from 0)
     *       (0, 2, 1, 3) dictates the starting indexes to which the matrix is added, 
     *       given the input <2,0> this means the indices of the 2nd (indicated by the "2") and the 4th (indicated by the "3") axes 
     *       are held fixed because they do not correspond to the entries <2,0>.
     */
    template <size_t r, size_t s, int Z = Rank>
    enable_if_t<Z == 4, Self&> addBlock(const MatrixX<Scalar>& M, size_t i, size_t j, size_t k, size_t l) {

        // Initialize series of arrays with 1 or 0 values, so that the correct tensor indices given by the template argument correspond to the matrix indices
        size_t ia[4] = {1, 0, 0, 0};
        size_t ja[4] = {0, 1, 0, 0};
        size_t ka[4] = {0, 0, 1, 0};
        size_t la[4] = {0, 0, 0, 1};

        for (size_t x = 0; x < M.rows(); x++) {
            for (size_t y = 0; y < M.cols(); y++) {

                size_t i_effective = i + x * ia[r] + y * ia[s];
                size_t j_effective = j + x * ja[r] + y * ja[s];
                size_t k_effective = k + x * ka[r] + y * ka[s];
                size_t l_effective = l + x * la[r] + y * la[s];

                this->operator()(i_effective, j_effective, k_effective, l_effective) += M(x, y);
            }
        }

        return (*this);
    }
};


}  // namespace GQCP
