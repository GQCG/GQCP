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


#include "Mathematical/Algorithm/Step.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"

#include <type_traits>


namespace GQCP {


/**
 *  An iteration step that calculates the matrix-vector products for all (new) guess vectors.
 */
class MatrixVectorProductCalculation:
    public Step<EigenproblemEnvironment> {

public:
    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Calculate the matrix-vector products for all the (new) guess vectors and add them to the environment.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(EigenproblemEnvironment& environment) override {

        const auto& V = environment.V;  // the subspace of guess vectors

        assert((V.transpose() * V).isApprox(MatrixX<double>::Identity(V.cols(), V.cols()), 1.0e-08));  // make sure that the subspace vectors are orthonormal


        auto& VA = environment.VA;  // VA = A * V (implicitly calculated through the matrix-vector product)
        const auto& matvec = environment.matrix_vector_product_function;

        // Check how many vectors there currently are in V and in VA: only calculate the expensive matrix-vector product for 'new' vectors.
        // If there is no difference, no matrix-vector products should be calculated.
        // If V is larger than VA:
        //      - VA should be expanded
        //      - we should calculate the matrix-vector products for the new vectors and place them into the columns of the new VA
        // If V is smaller than VA, which can only happen after a subspace collapse:
        //      - VA should shrink
        //      - we should calculate the matrix-vector products for all the vectors in the collapsed subspace
        const auto vectors_in_V = V.cols();
        const auto vectors_in_VA = VA.cols();
        const auto difference = vectors_in_V - vectors_in_VA;

        if (difference != 0) {
            VA.conservativeResize(Eigen::NoChange, VA.cols() + difference);  // accounts for both expansion and shrinking

            // Calculate the only the necessary matrix-vector products; find the start_index that accounts for both expansion and shrinking
            size_t start_index = 0;
            if ((difference > 0) && (vectors_in_VA > 0)) {
                start_index = vectors_in_VA - 1;  // -1 because of computers
            }

            for (size_t column_index = start_index; column_index < vectors_in_V; column_index++) {
                VA.col(column_index) = matvec(V.col(column_index));
            }
        }
    }
};


}  // namespace GQCP
