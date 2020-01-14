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


#include "Basis/TransformationMatrix.hpp"
#include "CISolver/CISolver.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "RDM/RDMCalculator.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  A class that is a wrapper around solving the eigenvalue problem for the molecular Hamiltonian using DOCI in the RHF basis
 */
class DOCIRHF {
private:
    Molecule molecule;  // the molecule that will be solved for
    std::string basis_set;  // the basisset that should be used

    bool is_solved = false;

    double energy_solution;
    double rhf_energy_solution;
    bool use_davidson;
    TransformationMatrix<double> T_total;  // total transformation from atomic orbital basis to the RHF orbitals

public:
    // CONSTRUCTORS

    /**
     *  @param xyz_filename         the file that contains the molecule specification (coordinates in angstrom)
     *  @param basis_set            the basisset that should be used
     */
    DOCIRHF(const std::string& xyz_filename, const std::string& basis_set, const bool use_davidson);


    /**
     *  @param molecule             the molecule that will be solved for
     *  @param basis_set            the basisset that should be used
     */
    DOCIRHF(const Molecule& molecule, const std::string& basis_set, const bool use_davidson);


    // PUBLIC METHODS

    /**
     *  Solve the eigenvalue problem for the molecular Hamiltonian in the doubly occupied space
     */
    void solve();

    /**
     *  @return DOCI energy
     */
    double energy() const;

    /**
     *  @return the RHF energy
     */
    double energy_rhf() const;

    /**
     *  @return the total transformation from atomic orbital basis to the RHF orbitals
     */
    const TransformationMatrix<double>& transformationMatrix() const;
};


}  // namespace QCMethod
}  // namespace GQCP
