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
#include "CISolver/CISolver.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "RDM/RDMCalculator.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  A class that diagonalizes the FCI Hamiltonian for a given molecule and performs Fukui and Dyson analysis.
 *  Fukui analysis revolves around the response of the 1RDM with respect to a change of N the number of electrons.
 *  Dyson analysis produces a Dyson orbital which is the overlap between the N-electron and N-1 electron wavefunction.
 *  The same orthonormal basis is used for both species (N and N-1 number of electrons): 
 *  RHF is used to determine an orthonormal basis for the species with an even number of electrons 
 */
class FukuiDysonAnalysis {
private:
    Molecule molecule;
    RSpinorBasis<double, GTOShell> sp_basis;
    SQHamiltonian<double> sq_hamiltonian;
    std::string basis_set;  // the basisset that should be used
    ProductFockSpace fock_space1 = ProductFockSpace(0, 0, 0);  // default
    ProductFockSpace fock_space2 = ProductFockSpace(0, 0, 0);  // default

    VectorX<double> dyson_coefficients;
    OneRDM<double> fukui_matrix;
    VectorX<double> fukui_naturals;
    SquareMatrix<double> fukui_vectors;

public:

    // CONSTRUCTORS
    /**
     *  @param molecule                 the molecule that will be solved for
     *  @param basis_set                the basisset that should be used
     *  @param use_diis                 flags if one wants to use the DIIS for the RHF solve as opposed to Plain RHF
     */
    FukuiDysonAnalysis(const Molecule& molecule, const std::string& basis_set, const bool use_diis);


    // GETTERS
    VectorX<double> get_dyson_coefficients() const { return this->dyson_coefficients; };
    OneRDM<double> get_fukui_matrix() const { return this->fukui_matrix; };
    VectorX<double> get_fukui_naturals() const { return this->fukui_naturals; };
    SquareMatrix<double> get_fukui_vectors() const { return this->fukui_vectors; };
    SquareMatrix<double> get_canonical_matrix() const { return this->sp_basis.transformationMatrix(); };
};


}  // namespace QCMethod
}  // namespace GQCP
