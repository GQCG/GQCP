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
#include "Mathematical/Representation/Tensor.hpp"

#include <map>
#include <numeric>
#include <vector>


namespace GQCP {


/**
 *  A slice of a rank-four tensor that only exists implicitly.
 * 
 *  If the full tensor is unnecessary to know, and only a certain slice of the tensor is of interest, this class implements operator() that can be used with the indices of the full tensor.
 */
template <typename _Scalar>
class ImplicitRankFourTensorSlice {
public:
    using Scalar = _Scalar;


private:
    std::vector<std::map<size_t, size_t>> indices_implicit_to_dense;  // an array of maps, mapping the implicit tensor indices to these of the dense representation

    Tensor<Scalar, 4> T;  // the dense representation of the slice


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Initialize an ImplicitRankFourTensorSlice's members.
     * 
     *  @param indices_implicit_to_dense                an array of maps, mapping the implicit tensor indices to these of the dense representation
     *  @param T                                        the dense representation of the slice
     */
    ImplicitRankFourTensorSlice(const std::vector<std::map<size_t, size_t>>& indices_implicit_to_dense, const Tensor<Scalar, 4>& T) :
        indices_implicit_to_dense {indices_implicit_to_dense},
        T {T} {

        const auto dimensions = T.dimensions();
        for (size_t axis_index = 0; axis_index < 4; axis_index++) {
            if (this->indices_implicit_to_dense[axis_index].size() != dimensions[axis_index]) {
                throw std::invalid_argument("ImplicitRankFourTensorSlice(const std::vector<std::map<size_t, size_t>>&, const Tensor<Scalar, 4>&): The given dense representation of the slice has an incompatible dimension for axis number " + std::to_string(axis_index) + ".");
            }
        }
    }


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Construct an ImplicitRankFourTensor from a dense tensor block and corresponding index ranges.
     * 
     *  @param axis1_start              the zero-based index of the implicit tensor at which the first axis should start
     *  @param axis1_end                the zero-based index of the implicit tensor at which the first axis should end (not included)
     *  @param axis2_start              the zero-based index of the implicit tensor at which the second axis should start
     *  @param axis2_end                the zero-based index of the implicit tensor at which the second axis should end (not included)
     *  @param axis3_start              the zero-based index of the implicit tensor at which the third axis should start
     *  @param axis3_end                the zero-based index of the implicit tensor at which the third axis should end (not included)
     *  @param axis4_start              the zero-based index of the implicit tensor at which the fourth axis should start
     *  @param axis4_end                the zero-based index of the implicit tensor at which the fourth axis should end (not included)
     *  @param T                        the dense representation of the block
     * 
     *  @return an implicit rank-four tensor slice
     */
    static ImplicitRankFourTensorSlice<Scalar> FromBlockRanges(const size_t axis1_start, const size_t axis1_end, const size_t axis2_start, const size_t axis2_end, const size_t axis3_start, const size_t axis3_end, const size_t axis4_start, const size_t axis4_end, const Tensor<Scalar, 4>& T) {

        // Convert the ranges into allowed indices of the implicit tensor, and then use another named constructor.
        std::vector<size_t> axis1_indices(axis1_end - axis1_start);
        std::iota(axis1_indices.begin(), axis1_indices.end(), axis1_start);

        std::vector<size_t> axis2_indices(axis2_end - axis2_start);
        std::iota(axis2_indices.begin(), axis2_indices.end(), axis2_start);

        std::vector<size_t> axis3_indices(axis3_end - axis3_start);
        std::iota(axis3_indices.begin(), axis3_indices.end(), axis3_start);

        std::vector<size_t> axis4_indices(axis4_end - axis4_start);
        std::iota(axis4_indices.begin(), axis4_indices.end(), axis4_start);

        return ImplicitRankFourTensorSlice<Scalar>::FromIndices(axis1_indices, axis2_indices, axis3_indices, axis4_indices, T);
    }


    /**
     *  A default constructor setting everything to zero.
     */
    ImplicitRankFourTensorSlice() :
        // Use a named constructor for the default initialization.
        ImplicitRankFourTensorSlice(ImplicitRankFourTensorSlice<Scalar>::ZeroFromIndices({}, {}, {}, {})) {}


    /**
     *  Create an implicit rank-four tensor slice through a dense representation of the slice.
     * 
     *  @param axes_indices                     an array where the i-th element represents the indices (in order) for the i-th tensor axis of the implicit tensor that the indices of the dense representation of the slice correspond to
     *  @param T                                the dense representation of the slice
     * 
     *  @return an implicit rank-four tensor slice
     */
    static ImplicitRankFourTensorSlice<Scalar> FromIndices(const std::vector<std::vector<size_t>>& axes_indices, const Tensor<Scalar, 4>& T) {

        // Loop over all indices of all axes, mapping them to the indices of the dense representation of the slice.
        std::vector<std::map<size_t, size_t>> indices_maps;
        for (size_t axis_index = 0; axis_index < 4; axis_index++) {
            size_t dense_index = 0;  // index in the dense representation of the slice

            std::map<size_t, size_t> axis_map {};  // the index mapping for the current axis index
            for (const auto& implicit_index : axes_indices[axis_index]) {
                axis_map[implicit_index] = dense_index;
                dense_index++;  // the dense indices are contiguous
            }
            indices_maps.push_back(axis_map);
        }

        return ImplicitRankFourTensorSlice<Scalar>(indices_maps, T);
    }


    /**
     *  Create an implicit rank-four tensor slice through a dense representation of the slice.
     * 
     *  @param axis1_indices                    the indices (in order) for the first tensor axis of the implicit tensor that the indices of the dense representation of the slice correspond to
     *  @param axis2_indices                    the indices (in order) for the second tensor axis of the implicit tensor that the indices of the dense representation of the slice correspond to
     *  @param axis3_indices                    the indices (in order) for the third tensor axis of the implicit tensor that the indices of the dense representation of the slice correspond to
     *  @param axis4_indices                    the indices (in order) for the fourth tensor axis of the implicit tensor that the indices of the dense representation of the slice correspond to
     *  @param T                                the dense representation of the slice
     * 
     *  @return an implicit rank-four tensor slice
     */
    static ImplicitRankFourTensorSlice<Scalar> FromIndices(const std::vector<size_t>& axis1_indices, const std::vector<size_t>& axis2_indices, const std::vector<size_t>& axis3_indices, const std::vector<size_t>& axis4_indices, const Tensor<Scalar, 4>& T) {

        // Wrap the separate indices for the axes into one array and use another named constructor.
        std::vector<std::vector<size_t>> axes_indices {axis1_indices, axis2_indices, axis3_indices, axis4_indices};

        return ImplicitRankFourTensorSlice<Scalar>::FromIndices(axes_indices, T);
    }


    /**
     *  Construct a zero ImplicitRankFourTensor from block index ranges.
     * 
     *  @param axis1_start              the zero-based index of the implicit tensor at which the first axis should start
     *  @param axis1_end                the zero-based index of the implicit tensor at which the first axis should end (not included)
     *  @param axis2_start              the zero-based index of the implicit tensor at which the second axis should start
     *  @param axis2_end                the zero-based index of the implicit tensor at which the second axis should end (not included)
     *  @param axis3_start              the zero-based index of the implicit tensor at which the third axis should start
     *  @param axis3_end                the zero-based index of the implicit tensor at which the third axis should end (not included)
     *  @param axis4_start              the zero-based index of the implicit tensor at which the fourth axis should start
     *  @param axis4_end                the zero-based index of the implicit tensor at which the fourth axis should end (not included)
     * 
     *  @return a zero implicit rank-four tensor slice
     */
    static ImplicitRankFourTensorSlice<Scalar> ZeroFromBlockRanges(const size_t axis1_start, const size_t axis1_end, const size_t axis2_start, const size_t axis2_end, const size_t axis3_start, const size_t axis3_end, const size_t axis4_start, const size_t axis4_end) {

        // Create the zero dense representation of the slice and then use another named constructor.
        const auto axis1_dimension = static_cast<long>(axis1_end - axis1_start);  // need static_cast for Tensor
        const auto axis2_dimension = static_cast<long>(axis2_end - axis2_start);
        const auto axis3_dimension = static_cast<long>(axis3_end - axis3_start);
        const auto axis4_dimension = static_cast<long>(axis4_end - axis4_start);

        Tensor<Scalar, 4> T {axis1_dimension, axis2_dimension, axis3_dimension, axis4_dimension};
        T.setZero();

        return ImplicitRankFourTensorSlice<Scalar>::FromBlockRanges(axis1_start, axis1_end, axis2_start, axis2_end, axis3_start, axis3_end, axis4_start, axis4_end, T);
    }


    /**
     *  Create a zero-initialized implicit rank-four tensor slice from allowed axes indices.
     * 
     *  @param axis1_indices                    the indices (in order) for the first tensor axis of the implicit tensor that the indices of the dense representation of the slice correspond to
     *  @param axis2_indices                    the indices (in order) for the second tensor axis of the implicit tensor that the indices of the dense representation of the slice correspond to
     *  @param axis3_indices                    the indices (in order) for the third tensor axis of the implicit tensor that the indices of the dense representation of the slice correspond to
     *  @param axis4_indices                    the indices (in order) for the fourth tensor axis of the implicit tensor that the indices of the dense representation of the slice correspond to
     * 
     *  @return a zero-initialized implicit rank-four tensor slice
     */
    static ImplicitRankFourTensorSlice<Scalar> ZeroFromIndices(const std::vector<size_t>& axis1_indices, const std::vector<size_t>& axis2_indices, const std::vector<size_t>& axis3_indices, const std::vector<size_t>& axis4_indices) {

        // Zero-initialize a tensor with the required dimensions and then use another named constructor.
        const auto axis1_dimension = static_cast<long>(axis1_indices.size());
        const auto axis2_dimension = static_cast<long>(axis2_indices.size());
        const auto axis3_dimension = static_cast<long>(axis3_indices.size());
        const auto axis4_dimension = static_cast<long>(axis4_indices.size());

        Tensor<Scalar, 4> T {axis1_dimension, axis2_dimension, axis3_dimension, axis4_dimension};
        T.setZero();

        return ImplicitRankFourTensorSlice<Scalar>::FromIndices(axis1_indices, axis2_indices, axis3_indices, axis4_indices, T);
    }


    /*
     *  OPERATORS
     */

    /**
     *  Access an element of this implicit tensor slice.
     * 
     *  @param index1           the first index of the implicit encapsulating tensor
     *  @param index2           the second index of the implicit encapsulating tensor
     *  @param index3           the third index of the implicit encapsulating tensor
     *  @param index4           the fourth index of the implicit encapsulating tensor
     * 
     *  @return a read-only element of the encapsulating tensor
     */
    Scalar operator()(const size_t index1, const size_t index2, const size_t index3, const size_t index4) const {

        // Map the implicit indices to the indices of the dense representation of this slice.
        const auto dense_index1 = this->denseIndexOf<0>(index1);
        const auto dense_index2 = this->denseIndexOf<1>(index2);
        const auto dense_index3 = this->denseIndexOf<2>(index3);
        const auto dense_index4 = this->denseIndexOf<3>(index4);

        return this->T(dense_index1, dense_index2, dense_index3, dense_index4);
    }

    /**
     *  Access an element of this implicit tensor slice.
     * 
     *  @param index1           the first index of the implicit encapsulating tensor
     *  @param index2           the second index of the implicit encapsulating tensor
     *  @param index3           the third index of the implicit encapsulating tensor
     *  @param index4           the fourth index of the implicit encapsulating tensor
     * 
     *  @return a writable element of the encapsulating tensor
     */
    Scalar& operator()(const size_t index1, const size_t index2, const size_t index3, const size_t index4) {

        // Map the implicit indices to the indices of the dense representation of this slice.
        const auto dense_index1 = this->denseIndexOf<0>(index1);
        const auto dense_index2 = this->denseIndexOf<1>(index2);
        const auto dense_index3 = this->denseIndexOf<2>(index3);
        const auto dense_index4 = this->denseIndexOf<3>(index4);

        return this->T(dense_index1, dense_index2, dense_index3, dense_index4);
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return this as a (column-major) matrix
     */
    MatrixX<Scalar> asMatrix() const { return this->T.pairWiseReduce(); }

    /**
     *  @return this as a tensor
     */
    const Tensor<Scalar, 4>& asTensor() const { return this->T; }

    /**
     *  Convert an implicit axis index to the axis index in the dense representation of this slice.
     * 
     *  @tparam Axis                the index of the axis (0,1,2,3)
     * 
     *  @param index                the implicit tensor index inside the given axis
     * 
     *  @return the index of the dense representation of this slice for the given axis
     */
    template <size_t Axis>
    size_t denseIndexOf(const size_t index) const { return this->indices_implicit_to_dense[Axis].at(index); }

    /**
     *  @return an array of maps, mapping the implicit tensor indices to these of the dense representation
     */
    const std::vector<std::map<size_t, size_t>>& indexMaps() const { return this->indices_implicit_to_dense; }
};


}  // namespace GQCP
