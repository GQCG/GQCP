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


#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/QCMethodProtocol.hpp"
#include "QCMethod/QCObjective.hpp"
#include "QCModel/HF/RHF.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  The restricted Hartree-Fock quantum chemical method
 * 
 *  @tparam _ExpansionScalar_       the type of scalar that is used for the expansion of the spatial orbitals in their underlying scalar basis
 */
template <typename _ExpansionScalar>
class RHF: public GQCP::QCMethodProtocol<QCModel::RHF<_ExpansionScalar>, QCMethod::RHF<_ExpansionScalar>> {
public:
    using ExpansionScalar = _ExpansionScalar;


private:
    size_t N_P;  // the number of electron pairs
    SQHamiltonian<ExpansionScalar> sq_hamiltonian;  // the Hamiltonian expressed in a scalar basis


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param N_P                  the number of electron pairs
     *  @param sq_hamiltonian       the Hamiltonian expressed in a scalar basis
     */
    RHF(const size_t N_P, const SQHamiltonian<ExpansionScalar>& sq_hamiltonian) :
        N_P (N_P),
        sq_hamiltonian (sq_hamiltonian)
    {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Optimize the electronic structure model: find the parameters that are the solutions to the quantum chemical method's objective
     * 
     *  @tparam QCObjective         the type of the objective
     *  @tparam Solver              the type of the solver
     * 
     *  @param objective            the objective that should be fulfilled in order to consider the model's parameters as 'optimal'
     *  @param solver               the solver that will try to optimize the parameters
     */
    template <typename QCObjective, typename Solver>
    QCStructure<QCModel::RHF<ExpansionScalar>> optimize(const QCObjective& objective, Solver& solver) {

        // In order to construct a valid QCModel::RHF, we need (1) the coefficient matrix, (2) the orbital energies
        // To make a QCStructure<QCModel::RHF<ExpansionScalar>>, we need the electronic energy, coefficient matrix, orbital energies and the number of electrons
        const auto C = solver.solve();  // the canonical coefficient matrix

        // Transform the Fock matrix to the orthonormal spinor basis: the orbital energies are on the diagonal
        const auto D = QCModel::RHF<ExpansionScalar>::calculateScalarBasis1RDM(C, 2*this->N_P);
        const auto F = QCModel::RHF<ExpansionScalar>::calculateScalarBasisFockMatrix(D, this->sq_hamiltonian);  // the converged Fock matrix (expressed in the scalar orbital basis)
        auto F_orthonormal = F;
        F_orthonormal.transform(C);  // now in the orthonormal spinor basis
        const auto orbital_energies = F_orthonormal.parameters().diagonal();

        // A PlainRHFSCFSolver only finds the ground state, so the QCStructure only needs to contain the parameters for one state
        const auto E_electronic = QCModel::RHF<ExpansionScalar>::calculateElectronicEnergy(D, this->sq_hamiltonian.core(), F);
        const QCModel::RHF<ExpansionScalar> rhf_parameters (N_P, orbital_energies, C);

        // Now that we have constructed an instance of the QCModel, we should check if the objective is fulfilled
        if (!objective.isSatisfiedWith(rhf_parameters)) {
            throw std::runtime_error("QCModel::RHF::optimize(const QCObjective&, Solver&): The solver did not produce a solution that fulfills the objective.");
        }
        return QCStructure<QCModel::RHF<ExpansionScalar>>({E_electronic}, {rhf_parameters});
    }
};


}  // namespace QCMethod
}  // namespace GQCP
