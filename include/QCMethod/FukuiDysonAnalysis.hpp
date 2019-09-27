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

#include "Basis/SingleParticleBasis.hpp"
#include "CISolver/CISolver.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "RDM/RDMCalculator.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  A class that solves the FCI Hamiltonian for a given molecule and performs Fukui and Dyson analysis
 *  restricted Hartree-Fock is initially performed for the singlet, 
 *  after which we solve in the FCI space and solve in the FCI space again after removing an electron
 */
class FukuiDysonAnalysis {
private:
    Molecule molecule;
    SingleParticleBasis<double, GTOShell> sp_basis;
    SQHamiltonian<double> sq_hamiltonian;
    std::string basis_set;  // the basisset that should be used
    ProductFockSpace fock_space1 = ProductFockSpace(0, 0, 0); // Default
    ProductFockSpace fock_space2 = ProductFockSpace(0, 0, 0); // Default

    VectorX<double> dyson_coefficients;
    OneRDM<double> fukui_matrix;
    OneRDM<double> fukui_naturals;
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
    OneRDM<double> get_fukui_naturals() const { return this->fukui_naturals; };
    SquareMatrix<double> get_fukui_vectors() const { return this->fukui_vectors; };
    SquareMatrix<double> get_canonical_matrix() const { return this->sp_basis.transformationMatrix(); };
};


}  // namespace QCMethod
}  // namespace GQCP
