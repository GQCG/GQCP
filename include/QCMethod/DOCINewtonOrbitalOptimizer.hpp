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
#ifndef GQCP_QCMETHOD_DOCINEWTONORBITALOPTIMIZER_HPP
#define GQCP_QCMETHOD_DOCINEWTONORBITALOPTIMIZER_HPP


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "RDM/RDMCalculator.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  A class that is a wrapper around solving the eigenvalue problem for the molecular Hamiltonian while optimizing the orbitals
 */
class DOCINewtonOrbitalOptimizer {
private:
    std::string xyz_filename;  // the file that contains the molecule specification (coordinates in angstrom)
    std::string basis_set;  // the basisset that should be used

    bool is_solved = false;
    bool use_davidson = false;
    bool localize = false;

    double energy_solution;
    SquareMatrix<double> T_total;  // total transformation from atomic orbital basis to the OO-DOCI orbitals

public:
    // CONSTRUCTORS

    /**
     *  @param xyz_filename         the file that contains the molecule specification (coordinates in angstrom)
     *  @param basis_set            the basisset that should be used
     *  @param use_davidson         indicate if one wants to use davidson to solve the eigenvalue problem (opposed to dense)
     *  @param localize             indicate if one wants to localize the orbitals before 
     */
    DOCINewtonOrbitalOptimizer(const std::string xyz_filename, const std::string basis_set, const bool use_davidson = false, const bool localize = false);


    // PUBLIC METHODS

    /**
     *  Solve the eigenvalue problem for the molecular Hamiltonian in the doubly occupied space
     */
    void solve();

    /**
     *  @return the newton orbital optimized ground state DOCI energy
     */
    double energy() const;

    /**
     *  @return the total transformation matrix to the OO-DOCI orbitals
     */
    SquareMatrix<double> transformationMatrix() const;
};


}  // namespace QCMethod
}  // namespace GQCP


#endif  // GQCP_QCMETHOD_DOCINEWTONORBITALOPTIMIZER_HPP
