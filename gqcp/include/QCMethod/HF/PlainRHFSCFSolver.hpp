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


#include "QCMethod/HF/BaseRHFSCFSolver.hpp"


namespace GQCP {


/**
 *  A class that implements a plain RHF SCF algorithm
 */
class PlainRHFSCFSolver : public BaseRHFSCFSolver {
private:
    /**
     *  Update the Fock matrix, i.e. calculate the Fock matrix to be used in the next iteration of the SCF procedure: the 'new' Fock matrix is just F = H_core + G
     *
     *  @param D_AO     the RHF density matrix in AO basis
     *
     *  @return the new Fock matrix (expressed in AO basis)
     */
    ScalarSQOneElectronOperator<double> calculateNewFockMatrix(const OneRDM<double>& D_AO) override;

public:
    // CONSTRUCTORS
    /**
     *  @param sq_hamiltonian                   the Hamiltonian expressed in an AO basis
     *  @param spinor_basis                     the spinor basis
     *  @param molecule                         the molecule used for the SCF calculation
     *  @param threshold                        the convergence treshold on the Frobenius norm on the AO density matrix
     *  @param maximum_number_of_iterations     the maximum number of iterations for the SCF procedure
     */
    PlainRHFSCFSolver(const SQHamiltonian<double>& sq_hamiltonian, const RSpinorBasis<double, GTOShell>& spinor_basis, const Molecule& molecule, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);
};


}  // namespace GQCP
