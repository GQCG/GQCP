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
#ifndef GQCP_EVALUATIONCONTAINER_HPP
#define GQCP_EVALUATIONCONTAINER_HPP


#include "math/SquareMatrix.hpp"

#include <Eigen/Sparse>


namespace GQCP {

/**
 *  Private class for Fock spaces
 */
template<class Container>
class EvaluationContainer {

    Container container;

    EvaluationContainer(size_t dimension) {
        container = Container::Zero(dimension, dimension);
    }

    void add(size_t i, size_t j, double value) {
        container(i, j) += value;
    }

    const Container& get_container() const {
        return container;
    }

    // Friend Classes
    friend class FockSpace;
    friend class SelectedFockSpace;
    friend class ProductFockSpace;
};


template<>
class EvaluationContainer<Eigen::SparseMatrix<double>> {

    Eigen::SparseMatrix<double> container;
    std::vector<Eigen::Triplet<double>> triplet_vector;

    EvaluationContainer(size_t dimension) {
        container = Eigen::SparseMatrix<double>(dimension, dimension);
    }

    void reserve(size_t n) {
        triplet_vector.reserve(n);
        container.reserve(n);
    }

    void add(size_t i, size_t j, double value) {
        this->triplet_vector.emplace_back(i, j, value);
    }

    void addToMatrix() {
        container.setFromTriplets(this->triplet_vector.begin(), this->triplet_vector.end());
        triplet_vector = {};
    }

    const Eigen::SparseMatrix<double>& get_container() const {
        return container;
    }

    const std::vector<Eigen::Triplet<double>> &get_triplets() const {
        return triplet_vector;
    };

    // Friend Classes
    friend class FockSpace;
    friend class SelectedFockSpace;
    friend class ProductFockSpace;
};


}  // namespace GQCP


#endif  // GQCP_EVALUATIONCONTAINER_HPP
