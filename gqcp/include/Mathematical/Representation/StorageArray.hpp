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


#include <vector>


namespace GQCP {


/**
 *  A simple array data structure.
 * 
 *  @tparam _Element        The type of element that this array stores.
 *  @tparam _Vectorizer     The type of the vectorizer that relates multiple tuple coordinates to a one-dimensional index.
 */
template <typename _Element, typename _Vectorizer>
// class StorageArray:
//     public VectorSpaceArithmetic<StorageArray<T, Vectorizer>> {
class StorageArray {
public:
    // The type of element that this array stores.
    using Element = _Element;

    // The type of the vectorizer that relates multiple tuple coordinates to a one-dimensional index.
    using Vectorizer = _Vectorizer;


private:
    // The one-dimensional representation of the elements of the array.
    std::vector<Element> m_elements;

    // The vectorizer that relates multiple tuple coordinates to a one-dimensional index.
    Vectorizer vectorizer;


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Create an array.
     * 
     *  @param elements         The one-dimensional representation of the elements of the array.
     *  @param vectorizer       The vectorizer that relates multiple tuple coordinates to a one-dimensional index.
     */
    StorageArray(const std::vector<Element>& elements, const Vectorizer& vectorizer) :
        m_elements {elements},
        vectorizer {vectorizer} {

        std::cout << "Created a vector of " << this->m_elements.size() << " elements whose data is at address " << this->m_elements.data() << std::endl;
    }


    /**
     *  Create an array with defaultly initialized elements.
     * 
     *  @param vectorizer       The vectorizer that relates multiple tuple coordinates to a one-dimensional index.
     */
    StorageArray(const Vectorizer& vectorizer) :
        StorageArray(std::vector<Element>(vectorizer.numberOfElements()), vectorizer) {}


    /*
     *  MARK: Element access
     */

    /**
     *  @param indices      A set of coordinates that accesses this array.
     * 
     *  @return The element that is stored at the given coordinate indices.
     */
    template <typename... Indices>
    Element operator()(const Indices&... indices) const {

        static_assert(sizeof...(indices) == Vectorizer::NumberOfAxes, "The number of indices must match the number of axes.");

        // Convert the indices pack to a vector so we can easily traverse.
        std::vector<size_t> indices_vector {static_cast<size_t>(indices)...};
        std::array<size_t, Vectorizer::NumberOfAxes> indices_array {};
        std::copy(indices_vector.begin(), indices_vector.end(), indices_array.begin());

        // Use the underlying vectorizer to produce the 1D offset, and access the 1D storage array accordingly.
        const auto vector_index = vectorizer.offset(indices_array);
        return this->m_elements[vector_index];
    }


    /**
     *  @param indices      A set of coordinates that accesses this array.
     * 
     *  @return The element that is stored at the given coordinate indices.
     */
    template <typename... Indices>
    Element& operator()(const Indices&... indices) {

        static_assert(sizeof...(indices) == Vectorizer::NumberOfAxes, "The number of indices must match the number of axes.");

        // Convert the indices pack to a vector so we can easily traverse.
        std::vector<size_t> indices_vector {static_cast<size_t>(indices)...};
        std::array<size_t, Vectorizer::NumberOfAxes> indices_array {};
        std::copy(indices_vector.begin(), indices_vector.end(), indices_array.begin());

        // Use the underlying vectorizer to produce the 1D offset, and access the 1D storage array accordingly.
        const auto vector_index = vectorizer.offset(indices_array);
        return this->m_elements[vector_index];
    }


    /**
     *  @return A read-only vector of all elements that this array stores.
     */
    const std::vector<Element>& elements() const { return this->m_elements; }


    /**
     *  @return A writable vector of all elements that this array stores.
     */
    std::vector<Element>& elements() { return this->m_elements; }


    // StorageArray<T, Vectorizer>& operator+=(const StorageArray<T, Vectorizer>& rhs) override {

    //     // TODO: Check if the vectorizers are actually equal.

    //     std::transform(this->elements.begin(), this->elements.end(),
    //                    rhs.elements.begin(), this->elements.begin(),
    //                    std::plus<T>());
    //     return (*this);
    // };


    // StorageArray<T, Vectorizer> operator-() const override {

    //     std::vector<T> result_elements;
    //     result_elements.reserve(this->elements.size());

    //     std::transform(this->elements.begin(), this->elements.end(),
    //                    std::back_inserter(result_elements), std::negate<T>());

    //     return StorageArray<T, Vectorizer> {result_elements, this->vectorizer};
    // }


    // Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> asEigenMatrix() const {
    //     using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    //     const auto rows = vectorizer.dimension(0);
    //     const auto cols = vectorizer.dimension(1);

    //     return Eigen::Map<const MatrixType>(this->elements.data(), rows, cols);
    // }
};


}  // namespace GQCP
