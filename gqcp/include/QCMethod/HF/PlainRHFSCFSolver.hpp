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
 *  A restricted Hartree-Fock self-consistent field (RHF SCF) solver.
 * 
 *  @tparam _ExpansionScalar        the type of scalar that is used to describe the expansion coefficients
 */
template <typename _ExpansionScalar>
class PlainRHFSCFSolver : public IterativeSolver<TransformationMatrix<_ExpansionScalar>, PlainRHFSCFSolver<_ExpansionScalar>> {
public:
    using ExpansionScalar = _ExpansionScalar;
    using Base = IterativeSolver<TransformationMatrix<_ExpansionScalar>, PlainRHFSCFSolver<_ExpansionScalar>>;

private:
    size_t N;  // the total number of electrons
    double threshold;  // the convergence threshold on the norm of the difference of consecutive density matrices
    SquareMatrix<ExpansionScalar> S;  // the overlap matrix of the spinor basis

    OneRDM<ExpansionScalar> D_previous;  // expressed in the scalar orbital basis
    OneRDM<ExpansionScalar> D_current;  // expressed in the scalar orbital basis
    SQHamiltonian<ExpansionScalar> sq_hamiltonian;  // expressed in the scalar orbital basis

public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  A general constructor from the properties
     * 
     *  @param C_initial                            the initial guess for the coefficient matrix
     *  @param N                                    the total number of electrons
     *  @param S                                    the overlap matrix of the spinor basis
     *  @param sq_hamiltonian                       the Hamiltonian expressed in that spinor basis
     *  @param threshold                            the convergence threshold on the norm of the difference of consecutive density matrices
     *  @param maximum_number_of_iterations         the maximum number of iterations the solver may perform
     */
    PlainRHFSCFSolver(const TransformationMatrix<ExpansionScalar>& C_initial, const size_t N, const SquareMatrix<ExpansionScalar>& S, const SQHamiltonian<ExpansionScalar>& sq_hamiltonian, const double threshold=1.0e-08, const size_t maximum_number_of_iterations=128) :
        Base(C_initial, maximum_number_of_iterations),
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
     *  NAMED CONSTRUCTORS
     */

    /**
     *  @param spinor_basis                         the spinor basis, containing the (non-orthogonal) underlying scalar orbitals
     *  @param sq_hamiltonian                       the Hamiltonian expressed in that spinor basis
     *  @param N                                    the total number of electrons
     *  @param threshold                            the convergence threshold on the norm of the difference of consecutive density matrices
     *  @param maximum_number_of_iterations         the maximum number of iterations the solver may perform
     * 
     *  @return an solver instance whose initial guess is obtained by diagonalizing the core Hamiltonian
     */
    static PlainRHFSCFSolver<ExpansionScalar> WithCoreGuess(const RSpinorBasis<ExpansionScalar, GTOShell>& spinor_basis, const SQHamiltonian<ExpansionScalar>& sq_hamiltonian, const size_t N, const double threshold=1.0e-08, const size_t maximum_number_of_iterations = 128) {

        const auto& H_core = sq_hamiltonian.core().parameters();
        const auto S = spinor_basis.overlap().parameters();

        // Obtain an initial guess for the density matrix in the scalar orbital basis by solving the generalized eigenvalue problem for H_core
        using MatrixType = Eigen::Matrix<ExpansionScalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver (H_core, S);
        TransformationMatrix<ExpansionScalar> C_initial = generalized_eigensolver.eigenvectors();

        return PlainRHFSCFSolver<ExpansionScalar>(C_initial, N, S, sq_hamiltonian, threshold, maximum_number_of_iterations);
    }


    /*
     *  PUBLIC ('OVERRIDDEN') METHODS
     */

    /**
     *  @return the converged Fock matrix (expressed in the scalar basis), i.e. calculated with converged density matrix, which is in turn calculated by the the coefficient matrix that is considered to be the solution
     */
    ScalarSQOneElectronOperator<ExpansionScalar> convergedFockMatrix() const {
        if (this->numberOfIterations() == 0) {
            throw std::runtime_error("PlainRHFSCFSolver::convergedFockMatrix(): You are trying to get the converged Fock matrix, but no iterations have been performed.");
        }

        return QCModel::RHF<ExpansionScalar>::calculateScalarBasisFockMatrix(this->D_current, this->sq_hamiltonian);
    }


    /**
     *  @return if the algorithm is considered to be converged
     */
    bool isConverged() {
        if (this->numberOfIterations() == 0) {
            return false;
        } else {
            return ((this->D_current - this->D_previous).norm() < this->threshold);
        } 
    }

    /**
     *  Diagonalize the Fock matrix (constructed from the current coefficient matrix) to find a next iteration of the coefficient matrix
     * 
     *  @return the next iterate
     */
    TransformationMatrix<ExpansionScalar> updateIterate() {
        const auto F = QCModel::RHF<ExpansionScalar>::calculateScalarBasisFockMatrix(this->D_current, this->sq_hamiltonian);  // the Fock matrix constructed from the current coefficient matrix

        // Solve the generalized eigenvalue problem for the Fock matrix to get a new iteration of the coefficient matrix
        using MatrixType = Eigen::Matrix<ExpansionScalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver (F.parameters(), this->S);
        const auto& C = generalized_eigensolver.eigenvectors();

        // Update the previous and current density matrices
        this->D_previous = this->D_current;
        this->D_current = QCModel::RHF<ExpansionScalar>::calculateScalarBasis1RDM(C, this->N);

        return C;
    }
};


}  // namespace GQCP
