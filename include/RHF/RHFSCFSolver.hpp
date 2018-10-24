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
 *  Base class for RHF SCF solvers. This class contains the solve()-method, which is the RHF SCF procedure.
 *
 *  Derived classes should implement the pure virtual function calculateNewFockMatrix.
 */
class RHFSCFSolver {
protected:
    const size_t maximum_number_of_iterations;
    const double threshold;
    bool is_converged = false;

    const GQCP::HamiltonianParameters ham_par;  // Hamiltonian parameters expressed in an AO basis
    const GQCP::Molecule molecule;

    GQCP::RHF solution;

    // PROTECTED METHODS
    /**
     *  Calculate a new Fock matrix (expressed in AO basis). This is the function that causes the different behaviour between derived classes.
     */
    virtual Eigen::MatrixXd calculateNewFockMatrix(const Eigen::MatrixXd& D_AO) = 0;

public:
    // CONSTRUCTORS
    /**
     *  Constructor based on given Hamiltonian parameters @param ham_par, @param molecule, @param maximum_number_of_iterations and @param SCF threshold
     */
    RHFSCFSolver(GQCP::HamiltonianParameters ham_par, GQCP::Molecule molecule, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);

    // GETTERS
    GQCP::RHF get_solution() const { return this->solution; }

    /**
     *  Solve the RHF SCF equations. This function internally uses the pure virtual calculateNewFockMatrix.
     */
    void solve();
};


}  // namespace GQCP


#endif /* RHFSCFSolver_hpp */
