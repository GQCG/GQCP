// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef RHFSCFSolver_hpp
#define RHFSCFSolver_hpp


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "RHF.hpp"
#include "Molecule.hpp"


namespace GQCP {


/**
 *  Base class for RHF SCF solvers. This class contains the solve()-method, which does the RHF SCF procedure.
 *
 *  Derived classes should implement the pure virtual function calculateNewFockMatrix().
 */
class RHFSCFSolver {
protected:
    size_t maximum_number_of_iterations;
    double threshold;
    bool is_converged = false;

    HamiltonianParameters ham_par;  // Hamiltonian parameters expressed in an AO basis
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
    virtual Eigen::MatrixXd calculateNewFockMatrix(const Eigen::MatrixXd& D_AO) = 0;

public:
    // CONSTRUCTORS
    /**
     *  @param ham_par                          the Hamiltonian parameters in AO basis
     *  @param molecule                         the molecule used for the SCF calculation
     *  @param threshold                        the convergence treshold on the Frobenius norm on the AO density matrix
     *  @param maximum_number_of_iterations     the maximum number of iterations for the SCF procedure
     */
    RHFSCFSolver(HamiltonianParameters ham_par, Molecule molecule, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);

    // GETTERS
    const RHF& get_solution() const { return this->solution; }

    /**
     *  Solve the RHF SCF equations
     */
    void solve();
};


}  // namespace GQCP


#endif /* RHFSCFSolver_hpp */
