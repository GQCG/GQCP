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


#include "Mathematical/Optimization/LinearEquation/LinearEquationEnvironment.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationSolver.hpp"
#include "Mathematical/Representation/Matrix.hpp"

#include <vector>


namespace GQCP {


/**
 *  An accelerator that uses a direct inversion of the iterative subspace (DIIS) on a subject to produce an accelerated subject.
 * 
 *  @tparam _Scalar             the scalar type that is used to represent an element of a DIIS error vector
 */
template <typename _Scalar>
class DIIS {
public:
    using Scalar = _Scalar;

public:
    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param subject              the subject
     * 
     *  @return the DIIS-accelerated subject
     */
    template <typename Subject>
    Subject accelerate(const std::vector<Subject>& subjects, const std::vector<VectorX<Scalar>>& errors) const {

        const auto diis_coefficients = this->calculateDIISCoefficients(errors);

        // Construct and return the DIIS-accelerated subject
        Subject accelerated_subject = diis_coefficients(0) * subjects.at(0);  // defaultly initializing may cause problems: the default constructor for a Matrix is a 0x0-matrix
        for (size_t i = 1; i < errors.size(); i++) {
            accelerated_subject += diis_coefficients(i) * subjects.at(i);
        }
        return accelerated_subject;
    }


    /**
     *  Find the linear combination of errors that minimizes the total error measure in the least squares sense (i.e. according to the DIIS algorithm).
     * 
     *  @return the coefficients that minimize the error measure
     */
    VectorX<Scalar> calculateDIISCoefficients(const std::vector<VectorX<Scalar>>& errors) const {

        const auto n = errors.size();

        // Initialize and calculate the augmented B matrix
        SquareMatrix<Scalar> B = -1 * SquareMatrix<Scalar>::Ones(n + 1, n + 1);  // +1 for the Lagrange multiplier
        B(n, n) = 0;
        for (size_t i = 0; i < n; i++) {
            const auto& error_i = errors[i];

            for (size_t j = 0; j < n; j++) {
                const auto& error_j = errors[j];
                B(i, j) = error_i.dot(error_j);
            }
        }

        // Initialize the RHS of the system of equations
        VectorX<Scalar> b = VectorX<Scalar>::Zero(n + 1);  // +1 for the multiplier
        b(n) = -1;                                         // the last entry of b is accessed through n: dimension of b is n+1 - 1 because of computers


        // Solve the DIIS linear equations [B x = b]
        auto environment = LinearEquationEnvironment<Scalar>(B, b);
        auto solver = LinearEquationSolver<Scalar>::HouseholderQR();
        solver.perform(environment);

        return environment.x;
    }
};


}  // namespace GQCP
