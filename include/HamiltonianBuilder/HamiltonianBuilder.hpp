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
#ifndef GQCP_HAMILTONIANBUILDER_HPP
#define GQCP_HAMILTONIANBUILDER_HPP


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "FockSpace/BaseFockSpace.hpp"

#include <memory>
#include <utility>



namespace GQCP {

/**
 *  Typedef for a type of function that handles where to 'put' a calculated value (e.g. matrix or vector).
 *
 *  @param I        index or address of an ONV
 *  @param J        index or address of ONV that couples with ONV I
 *  @param value    value related to the coupling
 */
using PassToMethod = std::function<void (size_t I, size_t J, double value)>;

/**
 *  A base class whose derived classes are able to construct matrix representations of the Hamiltonian in a Fock space
 *
 *  Derived classes should implement:
 *      - constructHamiltonian() which constructs the full Hamiltonian matrix in the given Fock space
 *      - matrixVectorProduct() which gives the result of the action of the Hamiltonian on a given coefficient vector
 *      - calculateDiagonal() which gives the diagonal of the Hamiltonian matrix
 */
class HamiltonianBuilder {
public:
    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor for the abstract base class
     */
    virtual ~HamiltonianBuilder() = 0;


    // PURE VIRTUAL GETTERS
    virtual const BaseFockSpace* get_fock_space() const = 0;


    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the Hamiltonian matrix
     */
    virtual Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters<double>& hamiltonian_parameters) const = 0;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *  @param x                            the vector upon which the Hamiltonian acts
     *  @param diagonal                     the diagonal of the Hamiltonian matrix
     *
     *  @return the action of the Hamiltonian on the coefficient vector
     */
    virtual Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters<double>& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) const = 0;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    virtual Eigen::VectorXd calculateDiagonal(const HamiltonianParameters<double>& hamiltonian_parameters) const = 0;
};


}  // namespace GQCP



#endif  // GQCP_HAMILTONIANBUILDER_HPP
