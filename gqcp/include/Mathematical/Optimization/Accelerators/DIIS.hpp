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


#include "Mathematical/Representation/Matrix.hpp"

#include <deque>


namespace GQCP {


/**
 *  An accelerator that uses a direct inversion of the iterative subspace (DIIS) on a subject to produce an accelerated subject.
 * 
 *  @tparam _Subject         the type whose instances should be accelerated
 */
template <typename _Subject, typename _Error, typename _ErrorScalar>
class DIIS {
public:
    using Subject = _Subject;
    using Error = _Error;
    using ErrorScalar = _ErrorScalar;


protected:
    size_t minimum_subspace_dimension;  // the minimum number of iterates that have to be in the subspace before enabling the DIIS acceleration
    size_t maximum_subspace_dimension;  // the maximum DIIS subspace dimension before the oldest Fock matrices get discarded (one at a time)


    std::deque<Subject> subjects;  // a deque that acts as a subspace of the available subjects
    std::deque<Error> errors;  // a deque that collects the corresponding error measures


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param minimum_subspace_dimension               the minimum number of subjects that have to be in the subspace before enabling the DIIS acceleration
     *  @param maximum_subspace_dimension               the maximum DIIS subspace dimension before the oldest subjects get discarded (one at a time)
     */
    DIIS(const size_t minimum_subspace_dimension=6, const size_t maximum_subspace_dimension=6) :
        minimum_subspace_dimension (minimum_subspace_dimension),
        maximum_subspace_dimension (maximum_subspace_dimension)
    {}


    /*
     *  PURE VIRTUAL PUBLIC METHODS
     */

    /**
     *  @param subject              the subject
     * 
     *  @return the error associated to the given subject
     */
    virtual Error calculateError(const Subject& subject) = 0;

    /**
     *  @param e1                   the first error
     *  @param e2                   the second error
     * 
     *  @return the scalar product between the two given errors
     */
    virtual ErrorScalar calculateErrorScalarProduct(const Error& e1, const Error& e2);


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Calculate an accelerated subject.
     * 
     *  @param subject              the subject
     * 
     *  @return an accelerated subject
     */
    Subject accelerate(const Subject& subject) {

        // Add the new subject to the subject subspace and perform an acceleration step if possible
        this->add_subject(subject);

        if (this->canAccelerate()) {
            const auto diis_coefficients = this->calculateDIISCoefficients();
            const auto accelerated_subject = this->calculateAcceleratedSubject(diis_coefficients);
            return accelerated_subject;
        }

        else {  // no acceleration is possible
            return subject;
        }
    }


    /**
     *  Add a subject to the subspace of subjects. Furthermore, calculate its associated error and manage the subspace of subjects if it becomes too large.
     * 
     *  @param subject              the subject
     */
    void add(const Subject& subject) {

        // Add the subject and its associated error to the subspace
        this->subjects.emplace(subject);
        const auto error = this->calculateError(subject);
        this->errors.emplace(error);

        // Remove the oldest subject and error if the subspace is large enough
        if (this->subjects.size() >= 2) {
            this->subjects.pop_front();
            this->errors.pop_front();
        }
    }


    /**
     *  @param coefficients             the coefficients that minimize the DIIS error (usually calculated from calculateDIISCoefficients())
     * 
     *  @return the DIIS-accelerated subject
     */
    Subject calculateAcceleratedSubject(const VectorX<ErrorScalar>& coefficients) {

        Subject accelerated_subject;  // the DIIS-accelerated subject
        for (size_t i = 0; i < this->errors.size(); i++) {
            accelerated_subject += coefficients(i) * this->subjects.at(i);  // a larger 'i' means a newer subject for .at(i)
        }
        return accelated_subject;
    }


    /**
     *  Minimize the error, according to the DIIS algorithm.
     * 
     *  @return the coefficients that minimize the error
     */
    VectorX<ErrorScalar> calculateDIISCoefficients() {

        const auto n = this->errors.size();

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
     *  @return if this accelerator can produce an accelerated subject
     */
    bool canAccelerate() const {
        if (this->subjects.size() >= this->minimum_subspace_dimension) {
            return true;
        } else {
            return false;
        }
    }
}


}  // namespace GQCP
