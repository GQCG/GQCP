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


#include <cstddef>
#include <numeric>


namespace GQCP {


/**
 *  An enumeration of possible axis ordering types.
 */
enum class Ordering {
    // Specifies that the first axis is contiguous.
    ColumnMajor,

    // Specifies that the last axis is contiguous.
    RowMajor
};


/**
 *  A class that takes care of mapping indices related to dense axes to vector indices.
 * 
 *  @tparam A   The number of axes supported by the vectorizer.
 */
template <size_t A>
class DenseVectorizer {
private:
    // The ordering of the axes.
    Ordering m_ordering;

    // The dimensions of each axis.
    std::array<size_t, A> m_dimensions;

    // The strides associated to each axis.
    std::array<size_t, A> m_strides;


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Initialize a `DenseVectorizer` from axis dimensions and an ordering type.
     * 
     *  @param dimensions       The dimensions of each axis.
     *  @param ordering         The ordering of the axes.
     */
    DenseVectorizer(const std::array<size_t, A>& dimensions, const Ordering ordering = Ordering::ColumnMajor) :
        m_ordering {ordering},
        m_dimensions {dimensions} {

        // Initialize the strides for every dimension.
        switch (this->m_ordering) {
        case Ordering::ColumnMajor: {  // the first axis is contiguous
            for (size_t k = 0; k < A; k++) {

                size_t stride {1};
                for (size_t l = 0; l < k; l++) {
                    stride *= this->m_dimensions[l];
                }

                this->m_strides[k] = stride;
            }
            break;
        }

        case Ordering::RowMajor: {  // the last axis is contiguous
            for (size_t k = 0; k < A; k++) {

                size_t stride {1};
                for (size_t l = k + 1; l < A; l++) {
                    stride *= this->m_dimensions[l];
                }

                this->m_strides[k] = stride;
            }
            break;
        }
        }
    }

    /**
     *  The default constructor.
     */
    DenseVectorizer() :
        DenseVectorizer(std::array<size_t, A> {}) {}


    /**
     *  MARK: Element access
     */

    /**
     *  Calculate the 1-dimensional offset number (also called: address, contiguous index) corresponding to the given indices that specify coordinates for each axis.
     * 
     *  @param indices      A list of indices, where each entry corresponds to the coordinate of each axis.
     * 
     *  @return The 1-dimensional offset number.
     */
    size_t offset(const std::array<size_t, A>& indices) const {

        size_t value {0};
        for (size_t k = 0; k < A; k++) {
            if (indices[k] > this->dimension(k)) {
                throw std::invalid_argument("DenseVectorizer::offset(const std::array<size_t, A>&): The index is not supported with the memory layout.");
            }

            value += this->stride(k) * indices[k];
        }

        return value;
    }


    /*
     *  MARK: General information
     */

    /**
     *  @param axis     An axis number.
     * 
     *  @return The dimension of the given axis.
     */
    size_t dimension(const size_t axis) const {
        if (axis >= A) {
            throw std::invalid_argument("DenseVectorizer::offset(const std::array<size_t, A>&): The given axis does not exist.");
        }

        return this->dimensions()[axis];
    }

    /**
     *  @return The dimensions of each axis.
     */
    const std::array<size_t, A>& dimensions() const { return this->m_dimensions; }

    // The number of axes for this vectorizer.
    static constexpr auto NumberOfAxes = A;

    /**
     *  @return The number of axes for this vectorizer.
     */
    size_t numberOfAxes() const { return A; }

    /**
     *  @return The number of elements that can be described by this vectorizer.
     */
    size_t numberOfElements() const {
        return std::accumulate(this->dimensions().begin(), this->dimensions().end(), 1, std::multiplies<size_t>());
    }


    /**
     *  @param axis     An axis number.
     * 
     *  @return The stride of the given axis.
     */
    size_t stride(const size_t axis) const {
        if (axis >= A) {
            std::cerr << "The given axis does not exist." << std::endl;
        }

        return this->strides()[axis];
    }


    /**
     *  @return The strides associated to each axis.
     */
    const std::array<size_t, A>& strides() const { return this->m_strides; }
};


/*
 *  MARK: Convenience aliases
 */

// A trivial vectorizer that offers tuple-based indexing for 0 arguments.
using ScalarVectorizer = DenseVectorizer<0>;

// A vectorizer that offers tuple-based indexing for 1 argument.
using VectorVectorizer = DenseVectorizer<1>;

// A vectorizer that offers tuple-based indexing for 2 arguments.
using MatrixVectorizer = DenseVectorizer<2>;

// A vectorizer that offers tuple-based indexing for a general number of arguments.
template <size_t A>
using TensorVectorizer = DenseVectorizer<A>;


}  // namespace GQCP
