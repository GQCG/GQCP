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
#ifndef GQCP_FROZENCOREINTERFACE_HPP
#define GQCP_FROZENCOREINTERFACE_HPP


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/HamiltonianBuilder.hpp"


namespace GQCP {


/**
 *  Base class containing general functions related to frozen core CI
 */
class FrozenCore : public HamiltonianBuilder {
protected:
    std::shared_ptr<HamiltonianBuilder> hamiltonian_builder;
    size_t X;

public:
    // CONSTRUCTORS
    FrozenCore(std::shared_ptr<HamiltonianBuilder>& hamiltonian_builder, size_t X);


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the Hamiltonian matrix
     */
    Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) const override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *  @param x                            the vector upon which the Hamiltonian acts
     *  @param diagonal                     the diagonal of the Hamiltonian matrix
     *
     *  @return the action of the Hamiltonian on the coefficient vector
     */
    Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) const override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    Eigen::VectorXd calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) const override;


    // PUBLIC METHODS
    /**
     *  @param hamiltonian_parameters              the Hamiltonian parameters in an orthonormal orbital basis
     *  @param X                                   amount of frozen core orbitals
     *
     *  @return new set of Hamiltonian parameters modified to fit the frozen core CI algorithm
     */
    HamiltonianParameters freezeHamiltonianParameters(const HamiltonianParameters& hamiltonian_parameters, size_t X) const;

    /**
     *  @param hamiltonian_parameters              the Hamiltonian parameters in an orthonormal orbital basis
     *  @param X                                   amount of frozen core orbitals
     *
     *  @return the diagonal correct value for the frozen core CI algorithm
     */
    double diagonalValue(const HamiltonianParameters& hamiltonian_parameters, size_t X) const;
};


}  // namespace GQCP


#endif //GQCP_FROZENCOREINTERFACE_HPP
