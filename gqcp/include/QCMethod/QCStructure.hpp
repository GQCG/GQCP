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
 *  A quantum chemical structure. It consists of the ground state (and possibly excited states) energy and the associated optimal quantum chemical parameters associated to a model.
 */
template <typename _QCModel>
class QCStructure {
public:
    using QCModel = _QCModel;


private:
    std::vector<double> energies;           // the ground (and possibly excited) state electronic energies
    std::vector<QCModel> model_parameters;  // the ground (and possibly excited) state model parameters


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param energies                 the ground (and possibly excited) state electronic energies
     *  @param model_parameters         the ground (and possibly excited) state model parameters
     */
    QCStructure(const std::vector<double>& energies, const std::vector<QCModel>& model_parameters) :
        energies {energies},
        model_parameters {model_parameters} {
        const auto n = energies.size();

        if (n < 1) {
            throw std::invalid_argument("QCStructure(const std::vector<double>&, const std::vector<QCModel>&): You have given an empty number of energies.");
        }

        if (n != model_parameters.size()) {
            throw std::invalid_argument("QCStructure(const std::vector<double>&, const std::vector<QCModel>&): The number of energies and sets of parameters do not match.");
        }
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param i            the index of the i-th excited state
     * 
     *  @return the electronic energy corresponding to the i-th excited state
     */
    double energy(const size_t i = 0) const { return this->energies[i]; }

    /**
     *  @return the ground state electronic energy for this quantum chemical structure
     */
    double groundStateEnergy() const { return this->energy(); }

    /**
     *  @return the ground state model parameters for this quantum chemical structure
     */
    const QCModel& groundStateParameters() const { return this->parameters(); }

    /**
     *  @param i            the index of the i-th excited state
     * 
     *  @return the parameters corresponding to the i-th excited state
     */
    const QCModel& parameters(const size_t i = 0) const { return this->model_parameters[i]; }
};


}  // namespace GQCP
