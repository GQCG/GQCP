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
#pragma once


#include "Mathematical/Optimization/IterativeSolver.hpp"

#include <deque>


namespace GQCP {


/**
 *  A base solver that accelerates its iterates using the direct inversion of the iterative subspace.
 * 
 *  @tparam _Iterate                the type of the iterate
 *  @tparam _Error                  the type of the error estimate
 *  @tparam _ErrorScalar            the scalar type of the scalar product between two errors
 */
template <typename _Iterate, typename _Error, typename _ErrorScalar>
class DIISSolver : public IterativeSolver<_Iterate> {
public:
    using Iterate = _Iterate;
    using Error = _Error;
    using ErrorScalar = _ErrorScalar;
    using Base = IterativeSolver<_Iterate>;


protected:
    size_t minimum_subspace_dimension;  // the minimum number of iterates that have to be in the subspace before enabling the DIIS acceleration
    size_t maximum_subspace_dimension;  // the maximum DIIS subspace dimension before the oldest Fock matrices get discarded (one at a time)

    std::deque<Iterate> iterate_deque;
    std::deque<Error> error_deque;


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param initial_guess                        the initial guess to the solver
     *  @param maximum_number_of_iterations         the maximum number of iterations the solver may perform   
     *  @param minimum_subspace_dimension       the minimum number of iterates that have to be in the subspace before enabling the DIIS acceleration
     *  @param maximum_subspace_dimension       the maximum DIIS subspace dimension before the oldest Fock matrices get discarded (one at a time)
     */
    DIISSolver(const Iterate& initial_guess, const size_t maximum_number_of_iterations=128, const size_t minimum_subspace_dimension=6, const size_t maximum_subspace_dimension=6) :
        Base(initial_guess, maximum_number_of_iterations),
        minimum_subspace_dimension (minimum_subspace_dimension),
        maximum_subspace_dimension (maximum_subspace_dimension)
    {}


    /*
     *  PUBLIC PURE VIRTUAL METHODS
     */

    /**
     *  @return a new iterate to be used in the next iteration if the DIIS acceleration has not been switched on yet
     */
    virtual Iterate regularNewIterate() = 0;

    /**
     *  @param iterate              the iterate
     * 
     *  @return the error associated to the given iterate
     */
    virtual Error calculateError(const Iterate& iterate) = 0;

    /**
     *  @param e1                   the first error
     *  @param e2                   the second error
     * 
     *  @return the scalar product between the two given errors
     */
    virtual ErrorScalar calculateErrorScalarProduct(const Error& e1, const Error& e2);


    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a new iterate according to the DIIS algorithm
     */
    virtual Iterate updateIterate() override {

        // Start by creating a regular new iterate and calculating the associated error
        Iterate iterate = this->regularNewIterate();
        Error error = this->calculateError(iterate);
        this->error_deque.emplace(error);


        // Apply the DIIS acceleration to produce a better iterate if the current subspace dimension is large enough
        if (this->subspaceDimension() >= this->minimum_subspace_dimension) {

            const auto diis_coefficients = this->calculateDIISCoefficients();
            Iterate diis_iterate;  // the DIIS-accelerated iterate
            for (size_t i = 0; i < this->subspaceDimension(); i++) {
                diis_iterate += diis_coefficients(i) * this->iterate_deque[i];
            }
            iterate = diis_iterate;  // the accelerated iterate should be used in the next iteration, and it should also be added to the current subspace
        }
        this->iterate_deque.emplace(iterate);


        // Discard the oldest iterate and errors if the subspace becomes too large
        if (this->subspaceDimension() > this->maximum_subspace_dimension) {
            this->iterate_deque.pop_front();
            this->error_deque.pop_front();
        }
        return iterate;
    }


    /*
     *  PUBLIC METHODS
     */

    VectorX<ErrorScalar> calculateDIISCoefficients() {

        const auto n = this->subspaceDimension();

        // Initialize and calculate the augmented B matrix
        SquareMatrix<ErrorScalar> B = -1 * SquareMatrix<ErrorScalar>::Ones(n+1,n+1);  // +1 for the multiplier
        B(n,n) = 0;
        for (size_t i = 0; i < n; i++) {
            const auto e_i = this->error_deque[i];

            for (size_t j = 0; j < n; j++) {
                const auto e_j = this->error_deque[j];
                B(i,j) = this->calculateErrorScalarProduct(e_i, e_j);
            }
        }

        // Initialize the RHS of the system of equations
        VectorX<ErrorScalar> b = VectorX<ErrorScalar>::Zero(n+1);  // +1 for the multiplier
        b(n) = -1;  // the last entry of b is accessed through n: dimension of b is n+1 - 1 because of computers


        // Solve the DIIS linear equations B y = b
        using MatrixType = Eigen::MatrixX<ErrorScalar, Eigen::Dynamic, Eigen::Dynamic;
        Eigen::HouseholderQR<MatrixType> linear_solver (B);
        return linear_solver.solve(b);
    }


    /**
     *  @return the current subspace dimension of the iterates
     */
    size_t subspaceDimension() const { return this->error_deque.size(); }
};


}  // namespace GQCP
