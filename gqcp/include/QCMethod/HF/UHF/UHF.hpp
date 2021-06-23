// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/UHF/UHFSCFEnvironment.hpp"
#include "QCMethod/QCStructure.hpp"
#include "QCModel/HF/UHF.hpp"

#include <type_traits>


namespace GQCP {
namespace QCMethod {


/**
 *  The unrestricted Hartree-Fock quantum chemical method.
 * 
 *  @tparam _Scalar             The type of scalar that is used for the expansion of the spatial orbitals in their underlying scalar basis: real or complex.
 */
template <typename _Scalar>
class UHF {

public:
    using Scalar = _Scalar;


public:
    /*
     *  PUBLIC METHODS
     */

    /**
     *  Optimize the UHF wave function model.
     * 
     *  @tparam Solver              The type of the solver.
     * 
     *  @param solver               The solver that will try to optimize the parameters.
     *  @param environment          The environment, which acts as a sort of calculation space for the solver.
     */
    template <typename Solver>
    QCStructure<QCModel::UHF<Scalar>, Scalar> optimize(Solver& solver, UHFSCFEnvironment<Scalar>& environment) const {

        // The UHF method's responsibility is to try to optimize the parameters of its method, given a solver and associated environment.
        solver.perform(environment);

        // To make a QCStructure<QCModel::UHF<Scalar>>, we need the electronic energy, coefficient matrices, orbital energies and the numbers of electrons.
        // Furthermore, the current UHF SCF solvers only find the ground state wave function parameters, so the QCStructure only needs to contain the parameters for one state.
        const auto& E_electronic = environment.electronic_energies.back();
        const auto& C = environment.coefficient_matrices.back();
        const auto& orbital_energies = environment.orbital_energies.back();
        const auto& N = environment.N;

        const QCModel::UHF<Scalar> uhf_parameters {N.alpha(), N.beta(), orbital_energies.alpha(), orbital_energies.beta(), C};

        return QCStructure<QCModel::UHF<Scalar>, Scalar>({E_electronic}, {uhf_parameters});
    }
};


}  // namespace QCMethod
}  // namespace GQCP
