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


namespace GQCP {


/**
 *  A step that adds projected correction vectors to the subspace (if their norm is large enough) and collapses the subspace if it becomes too large.
 */
class SubspaceUpdate:
    public Step<EigenproblemEnvironment> {


private:
    size_t maximum_subspace_dimension;
    double threshold;  // the threshold on the norm used for determining if a new projected correction vector should be added to the subspace


public:
    /*
     * CONSTRUCTORS
     */

    /**
     *  @param maximum_subspace_dimension           the maximum dimension of the subspace before collapsing
     *  @param threshold                            the threshold on the norm used for determining if a new projected correction vector should be added to the subspace
     */
    SubspaceUpdate(const size_t maximum_subspace_dimension = 15, const double threshold = 1.0e-03) :
        maximum_subspace_dimension {maximum_subspace_dimension},
        threshold {threshold} {}


    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Add projected correction vectors to the subspace (if their norm is large enough) and collapse the subspace becomes too large. The new subspace vectors (after a collapse) are linear combinations of current subspace vectors, with coefficients found in the lowest eigenvectors of the subspace matrix.";
    }


    /**
     *  Add projected correction vectors to the subspace (if their norm is large enough) and collapse the subspace becomes too large. The new subspace vectors (after a collapse) are linear combinations of current subspace vectors, with coefficients found in the lowest eigenvectors of the subspace matrix.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(EigenproblemEnvironment& environment) override {

        auto& V = environment.V;
        const auto& Delta = environment.Delta;

        // If the subspace will potentially become too large, collapse it in advance.
        const auto current_subspace_dimension = V.cols();
        if (current_subspace_dimension + Delta.cols() > this->maximum_subspace_dimension) {
            V = environment.X;
        }

        // Update the current subspace V with new vectors: add the normalized orthogonal projection of the correction vectors if their norm is large enough.
        // Note that we can't add more than one vector simultaneously, as the inclusion of one vector changes the subspace, which in turn changes its orthogonal complement.
        for (size_t column_index = 0; column_index < Delta.cols(); column_index++) {
            VectorX<double> v = Delta.col(column_index) - V * (V.transpose() * Delta.col(column_index));  // project the correction vector on the orthogonal complement of V
            const double norm = v.norm();
            v.normalize();

            if (norm > this->threshold) {
                V.conservativeResize(Eigen::NoChange, V.cols() + 1);  // the number of rows doesn't change
                V.col(V.cols() - 1) = v;                              // add the new vector to the last column
            }
        }
    }
};


}  // namespace GQCP
