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


#include "QCMethod/CC/CCSDEnvironment.hpp"
#include "QCMethod/QCStructure.hpp"
#include "QCModel/CC/CCD.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  The CCD quantum chemical method.
 * 
 *  @tparam _Scalar                 the scalar type used to represent the T2-amplitudes
 */
template <typename _Scalar>
class CCD {
public:
    using Scalar = _Scalar;

public:
    /*
     *  PUBLIC METHODS
     */

    /**
     *  Optimize the CCD wave function model.
     * 
     *  @tparam Solver              the type of the solver
     * 
     *  @param solver               the solver that will try to optimize the parameters
     *  @param environment          the environment, which acts as a sort of calculation space for the solver
     */
    template <typename Solver>
    QCStructure<GQCP::QCModel::CCD<Scalar>, Scalar> optimize(Solver& solver, CCSDEnvironment<Scalar>& environment) const {

        // The CCD method's responsibility is to try to optimize the parameters of its method, given a solver and associated environment.
        solver.perform(environment);

        // To make a QCStructure, we need the electronic (correlation) energy and the T2-amplitudes.
        // Furthermore, the solvers only find the ground state wave function parameters, so the QCStructure only needs to contain the parameters for one state.
        const auto& T2 = environment.t2_amplitudes.back();

        const auto E_electronic_correlation = environment.correlation_energies.back();
        const QCModel::CCD<Scalar> ccd_parameters {T2};

        return QCStructure<GQCP::QCModel::CCD<Scalar>, Scalar>({E_electronic_correlation}, {ccd_parameters});
    }
};


}  // namespace QCMethod
}  // namespace GQCP
