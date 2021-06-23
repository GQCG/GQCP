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


#include <stdexcept>

#include <vector>


namespace GQCP {


/**
 *  A quantum chemical structure. It encapsulates (energy, optimal parameter)-pairs for the ground state, and possibly excited states.
 * 
 *  @tparam _QCModel                The type of the optimal model parameters.
 *  @tparam _Scalar                 The scalar type of the energy: real or complex.
 */
template <typename _QCModel, typename _Scalar = double>
class QCStructure {
public:
    // The type of the optimal model parameters.
    using QCModel = _QCModel;

    // The scalar type of the energy: real or complex.
    using Scalar = _Scalar;


private:
    // The electronic energies for the ground state and possibly for the excited states.
    std::vector<Scalar> energies;

    // The optimal model parameters for the ground state and possibly for the excited states.
    std::vector<QCModel> model_parameters;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param energies                 The electronic energies for the ground state and possibly for the excited states.
     *  @param model_parameters         The optimal model parameters for the ground state and possibly for the excited states.
     */
    QCStructure(const std::vector<Scalar>& energies, const std::vector<QCModel>& model_parameters) :
        energies {energies},
        model_parameters {model_parameters} {

        const auto n = energies.size();

        if (n < 1) {
            throw std::invalid_argument("QCStructure(const std::vector<Scalar>&, const std::vector<QCModel>&): You have given an empty number of energies.");
        }

        if (n != model_parameters.size()) {
            throw std::invalid_argument("QCStructure(const std::vector<Scalar>&, const std::vector<QCModel>&): The number of energies and sets of parameters do not match.");
        }
    }


    /*
     *  MARK: Energy
     */

    /**
     *  @param i            The index of an excited state.
     * 
     *  @return The electronic energy corresponding to the i-th excited state.
     */
    Scalar energy(const size_t i = 0) const { return this->energies[i]; }

    /**
     *  @return The ground state electronic energy for this quantum chemical structure.
     */
    Scalar groundStateEnergy() const { return this->energy(); }


    /*
     *  MARK: Parameters
     */

    /**
     *  @param i            The index of an excited state.
     * 
     *  @return The optimal model parameters corresponding to the i-th excited state.
     */
    const QCModel& parameters(const size_t i = 0) const { return this->model_parameters[i]; }

    /**
     *  @return The ground state model parameters for this quantum chemical structure.
     */
    const QCModel& groundStateParameters() const { return this->parameters(); }
};


}  // namespace GQCP
