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


#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Optimization/IterativeSolver.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "QCModel/HF/RHF.hpp"

#include <Eigen/Dense>


namespace GQCP {

/**
 *  A restricted Hartree-Fock self-consistent field (RHF SCF) solver that uses a DIIS accelerator on the Fock matrix.
 * 
 *  This is not a regular DIIS solver, as it does not construct a linear combination of the iterates (i.e. the coefficient matrices): it uses a linear combination of previous Fock matrices, from which the next iterate is calculated.
 * 
 *  @tparam _ExpansionScalar        the type of scalar that is used to describe the expansion coefficients
 */
template <typename _ExpansionScalar>
class DIISRHFSCFSolver : public DIISSolver<TransformationMatrix<_ExpansionScalar>, SquareMatrix<_ExpansionScalar>, _ExpansionScalar> {
public:
    using ExpansionScalar = _ExpansionScalar;
    using Base = DIISSolver<TransformationMatrix<_ExpansionScalar>, SquareMatrix<_ExpansionScalar>, _ExpansionScalar>;


private:
    size_t N;  // the total number of electrons
    double threshold;  // the convergence threshold on the norm of the difference of consecutive density matrices
    SquareMatrix<ExpansionScalar> S;  // the overlap matrix of the spinor basis

    OneRDM<ExpansionScalar> D_previous;  // expressed in the scalar orbital basis
    OneRDM<ExpansionScalar> D_current;  // expressed in the scalar orbital basis
    SQHamiltonian<ExpansionScalar> sq_hamiltonian;  // expressed in the scalar orbital basis


public:

    /*
     * CONSTRUCTORS
     */

    /**
     *  A general constructor from the properties
     * 
     *  @param C_initial                            the initial guess for the coefficient matrix
     *  @param N                                    the total number of electrons
     *  @param S                                    the overlap matrix of the spinor basis
     *  @param sq_hamiltonian                       the Hamiltonian expressed in that spinor basis
     *  @param maximum_number_of_iterations         the maximum number of iterations the solver may perform   
     *  @param minimum_subspace_dimension           the minimum number of iterates that have to be in the subspace before enabling the DIIS acceleration
     *  @param threshold                            the convergence threshold on the norm of the difference of consecutive density matrices
     *  @param maximum_number_of_iterations         the maximum number of iterations the solver may perform
     */
    DIISRHFSCFSolver(const TransformationMatrix<ExpansionScalar>& C_initial, const size_t N, const SquareMatrix<ExpansionScalar>& S, const SQHamiltonian<ExpansionScalar>& sq_hamiltonian, const size_t minimum_subspace_dimension=6, const size_t maximum_subspace_dimension=6, const double threshold=1.0e-08, const size_t maximum_number_of_iterations=128) :
        Base(C_initial, maximum_number_of_iterations, minimum_subspace_dimension, maximum_subspace_dimension),
        threshold (threshold),
        N (N),
        S (S),
        sq_hamiltonian (sq_hamiltonian)
    {
        // Given the initial coefficient matrix, we can calculate the initial density matrix
        this->D_previous = QCModel::RHF<ExpansionScalar>::calculateScalarBasis1RDM(this->iterate, this->N);
        this->D_current = this->D_previous;  // at the start of the algorithm, they are equal
    }



    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a new iterate to be used in the next iteration if the DIIS acceleration has not been switched on yet
     */
    Iterate regularNewIterate() override {

        const auto F = QCModel::RHF<ExpansionScalar>::calculateScalarBasisFockMatrix(this->D_current, this->sq_hamiltonian);  // the Fock matrix constructed from the current coefficient matrix

    }


};


}  // namespace GQCP
