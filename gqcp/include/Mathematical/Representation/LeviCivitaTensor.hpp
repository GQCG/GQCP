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


#include "Mathematical/Representation/Tensor.hpp"

#include <unsupported/Eigen/CXX11/TensorSymmetry>


namespace GQCP {


/**
 *  The antisymmetric rank-three Levi-Civita tensor.
 * 
 *  @tparam _Scalar         The scalar type of one of the elements of the tensor. 
 */
template <typename _Scalar>
class LeviCivitaTensor {
public:
    // The scalar type of one of the elements of the tensor.
    using Scalar = _Scalar;

private:
    // The values of the Levi-Civita tensor.
    Tensor<Scalar, 3> epsilon;

public:
    /*
     *  MARK: Constructors
     */

    /**
     *  The default constructor.
     */
    LeviCivitaTensor() {

        // Initialize a zero tensor.
        this->epsilon = Tensor<Scalar, 3>(3, 3, 3);
        this->epsilon.setZero();


        // Construct the antisymmetries of the Levi-Civita tensor and set its elements.
        Eigen::SGroup<Eigen::AntiSymmetry<0, 1>, Eigen::AntiSymmetry<1, 2>> symmetry;
        symmetry(this->epsilon, 0, 1, 2) = 1.0;
    }


    /*
     *  MARK: Element access
     */

    /**
     *  Access one of the elements of the Levi-Civita tensor.
     *
     *  @param i            The first accessor index.
     *  @param j            The second accessor index.
     *  @param k            The third accessor index.
     * 
     *  @return The element epsilon(i,j,k) of the Levi-Civita tensor.
     */
    Scalar operator()(const size_t i, const size_t j, const size_t k) const { return this->epsilon(i, j, k); }

    /**
     *  Find the index k such that the Levi-Civita tensor epsilon(i,j,k) (or any permutations thereof) does not vanish.
     * 
     *  @param i            An index of the Levi-Civita tensor.
     *  @param j            Another index of the Levi-Civita tensor.
     * 
     *  @return The index k such that epsilon(i,j,k) does not vanish.
     */
    size_t nonZeroIndex(const size_t i, const size_t j) const {

        if (i == j) {
            throw std::invalid_argument("LeviCivitaTensor::nonZeroIndex(const size_t, const size_t): The given indices cannot be equal.");
        }

        return 3 - (i + j);
        
    }
};


}  // namespace GQCP
