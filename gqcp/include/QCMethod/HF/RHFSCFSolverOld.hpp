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
#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/RHF.hpp"


namespace GQCP {


/**
 *  Base class for RHF SCF solvers. This class contains the solve()-method, which does the RHF SCF procedure.
 *
 *  Derived classes should implement the pure virtual function calculateNewFockMatrix().
 */
class RHFSCFSolverOld {
protected:
    size_t maximum_number_of_iterations;
    double threshold;
    bool is_converged = false;

    RSpinorBasis<double, GTOShell> spinor_basis;  // the spinor basis (AOs)
    SQHamiltonian<double> sq_hamiltonian;  // the Hamiltonian expressed in an AO basis
    Molecule molecule;

    RHF solution;


    // PROTECTED METHODS

    /**
     *  Update the Fock matrix, i.e. calculate the Fock matrix to be used in the next iteration of the SCF procedure
     *
     *  @param D_AO     the RHF density matrix in AO basis
     *
     *  @return the new Fock matrix (expressed in AO basis)
     */
    virtual ScalarSQOneElectronOperator<double> calculateNewFockMatrix(const OneRDM<double>& D_AO) = 0;


public:
    // CONSTRUCTORS

    /**
     *  @param sq_hamiltonian                   the Hamiltonian expressed in an AO basis
     *  @param spinor_basis                     the spinor basis
     *  @param molecule                         the molecule used for the SCF calculation
     *  @param threshold                        the convergence treshold on the Frobenius norm on the AO density matrix
     *  @param maximum_number_of_iterations     the maximum number of iterations for the SCF procedure
     */
    RHFSCFSolverOld(const SQHamiltonian<double>& sq_hamiltonian, const RSpinorBasis<double, GTOShell>& spinor_basis, const Molecule& molecule, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);


    // DESTRUCTOR
    virtual ~RHFSCFSolverOld() = default;


    // GETTERS
    const RHF& get_solution() const { return this->solution; }


    // PUBLIC METHODS

    /**
     *  Solve the RHF SCF equations, obtaining an initial guess by solving the generalized eigenvalue problem for H_core
     */
    void solve();

    /**
     *  Solve the RHF SCF equations using an initial guess
     * 
     *  @param C_initial            the initial guess for the canonical RHF coefficient matrix
     */
    void solve(const TransformationMatrix<double>& C_initial);
};


}  // namespace GQCP
