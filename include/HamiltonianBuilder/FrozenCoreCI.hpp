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
#ifndef GQCP_FROZENCORE_HPP
#define GQCP_FROZENCORE_HPP


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/HamiltonianBuilder.hpp"


namespace GQCP {


/**
 *  (base) Class implementing general functions related to frozen core CI
 */
class FrozenCoreCI : public HamiltonianBuilder {
protected:
    size_t X;  // number of frozen orbitals/electrons
    std::shared_ptr<HamiltonianBuilder> active_hamiltonian_builder;  // non-frozen core Hamiltonian builder performing the HamiltonianBuilder interface in the active space with the frozen Hamiltonian parameters

public:
    // CONSTRUCTORS
    /**
     *  @param hamiltonian_builder           shared pointer to active (non-frozen core) Hamiltonian builder
     *  @param X                             the number of frozen orbitals
     */
    FrozenCoreCI(std::shared_ptr<HamiltonianBuilder> hamiltonian_builder, size_t X);


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the frozen core Hamiltonian matrix
     */
    Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) const override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *  @param x                            the vector upon which the Hamiltonian acts
     *  @param diagonal                     the diagonal of the Hamiltonian matrix
     *
     *  @return the action of the frozen core Hamiltonian on the coefficient vector
     */
    Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) const override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the frozen core Hamiltonian
     */
    Eigen::VectorXd calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) const override;


    // PUBLIC METHODS
    /**
     *  @param hamiltonian_parameters              the Hamiltonian parameters in an orthonormal orbital basis
     *  @param X                                   the number of frozen orbitals
     *
     *  @return a set of 'frozen' Hamiltonian parameters which cover two-electron integral evaluations from the active and inactive orbitals
     *  (see https://drive.google.com/file/d/1Fnhv2XyNO9Xw9YDoJOXU21_6_x2llntI/view?usp=sharing)
     */
    HamiltonianParameters freezeHamiltonianParameters(const HamiltonianParameters& hamiltonian_parameters, size_t X) const;

    /**
     *  @param hamiltonian_parameters              the Hamiltonian parameters in an orthonormal orbital basis
     *  @param X                                   the number of frozen orbitals
     *
     *  @return the diagonal from strictly evaluating the frozen orbitals in the Fock space
     */
    Eigen::VectorXd calculateFrozenCoreDiagonal(const HamiltonianParameters& hamiltonian_parameters, size_t X) const;
};


}  // namespace GQCP


#endif  // GQCP_FROZENCORE_HPP
