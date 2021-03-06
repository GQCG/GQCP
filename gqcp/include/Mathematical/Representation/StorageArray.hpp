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


#include "Mathematical/Representation/DenseVectorizer.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Utilities/type_traits.hpp"

#include <vector>


namespace GQCP {


/**
 *  A simple array data structure.
 * 
 *  @tparam _Element        The type of element that this array stores.
 *  @tparam _Vectorizer     The type of the vectorizer that relates multiple tuple coordinates to a one-dimensional index.
 */
template <typename _Element, typename _Vectorizer>
class StorageArray {
public:
    // The type of element that this array stores.
    using Element = _Element;

    // The type of the vectorizer that relates multiple tuple coordinates to a one-dimensional index.
    using Vectorizer = _Vectorizer;

    // The type of this;
    using Self = StorageArray<Element, Vectorizer>;

    // The number of axes for the underlying vectorizer.
    static constexpr auto NumberOfAxes = Vectorizer::NumberOfAxes;


private:
    // The one-dimensional representation of the elements of the array.
    std::vector<Element> m_elements;

    // The vectorizer that relates multiple tuple coordinates to a one-dimensional index.
    Vectorizer m_vectorizer;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create an array.
     * 
     *  @param elements         The one-dimensional representation of the elements of the array.
     *  @param vectorizer       The vectorizer that relates multiple tuple coordinates to a one-dimensional index.
     */
    StorageArray(const std::vector<Element>& elements, const Vectorizer& vectorizer) :
        m_elements {elements},
        m_vectorizer {vectorizer} {}


    /**
     *  Create an array with equal elements.
     * 
     *  @param elements         The one-dimensional representation of the elements of the array.
     *  @param vectorizer       The vectorizer that relates multiple tuple coordinates to a one-dimensional index.
     */
    StorageArray(const Element& element, const Vectorizer& vectorizer) :
        StorageArray(std::vector<Element>(vectorizer.numberOfElements(), element), vectorizer) {}


    /**
     *  Create an array with defaultly initialized elements.
     * 
     *  @param vectorizer       The vectorizer that relates multiple tuple coordinates to a one-dimensional index.
     */
    StorageArray(const Vectorizer& vectorizer) :
        StorageArray(std::vector<Element>(vectorizer.numberOfElements()), vectorizer) {}


    /**
     *  Create an array from a `std::array` and a vectorizer.
     * 
     *  @param elements         The one-dimensional representation of the elements of the array.
     *  @param vectorizer       The vectorizer that relates multiple tuple coordinates to a one-dimensional index.
     * 
     *  @tparam N               The number of elements in the `std::array`.
     */
    template <size_t N>
    StorageArray(const std::array<Element, N>& elements, const Vectorizer& vectorizer) :
        StorageArray(std::vector<Element>(elements.begin(), elements.end()), vectorizer) {}  // Convert the std::array into a std::vector.


    /*
     *  MARK: Conversions
     */

    /**
     *  A conversion operator between a scalar array and the underlying scalar.
     * 
     *  For the enable_if-part, check https://stackoverflow.com/a/18101382.
     */
    template <typename Z = Vectorizer, typename = typename std::enable_if<std::is_same<Z, ScalarVectorizer>::value>::type>
    operator Element() const {
        return (*this)();  // the call operator with no indices
    }


    /*
     *  MARK: Element access
     */

    /**
     *  @param indices      A set of coordinates that accesses this array.
     * 
     *  @return The element that is stored at the given coordinate indices.
     */
    template <typename... Indices>
    const Element& operator()(const Indices&... indices) const {

        static_assert(sizeof...(indices) == NumberOfAxes, "The number of indices must match the number of axes.");

        // Convert the indices pack to a vector so we can easily traverse.
        std::vector<size_t> indices_vector {static_cast<size_t>(indices)...};
        std::array<size_t, NumberOfAxes> indices_array {};
        std::copy(indices_vector.begin(), indices_vector.end(), indices_array.begin());

        // Use the underlying vectorizer to produce the 1D offset, and access the 1D storage array accordingly.
        const auto vector_index = this->vectorizer().offset(indices_array);
        return this->elements()[vector_index];
    }


    /**
     *  @param indices      A set of coordinates that accesses this array.
     * 
     *  @return The element that is stored at the given coordinate indices.
     */
    template <typename... Indices>
    Element& operator()(const Indices&... indices) { return const_cast<Element&>(const_cast<const Self*>(this)->operator()(indices...)); }


    /**
     *  @return A read-only vector of all elements that this array stores.
     */
    const std::vector<Element>& elements() const { return this->m_elements; }


    /**
     *  @return A writable vector of all elements that this array stores.
     */
    std::vector<Element>& elements() { return this->m_elements; }


    /*
     *  MARK: General information
     */

    /**
     *  @return The vectorizer that relates multiple tuple coordinates to a one-dimensional index.
     */
    const Vectorizer& vectorizer() const { return this->m_vectorizer; }


    /*
     *  MARK: Conversions
     */

    /**
     *  Convert this `StorageArray` to a `Vector`, if possible.
     * 
     *  @note This method is only available for `StorageArrays` that admit vector-like storage.
     */
    template <typename Z = Vectorizer>
    enable_if_t<std::is_same<Z, VectorVectorizer>::value, VectorX<Element>> asVector() {  // Marking this method `const` causes compile-time errors related to Eigen.

        return Eigen::Map<Eigen::Matrix<Element, Eigen::Dynamic, 1>> {this->elements().data(), static_cast<long>(this->vectorizer().dimension(0))};  // The dimension of the '0'-th axis is the number of rows, i.e. the dimension of the columns.
    }
};


}  // namespace GQCP
