// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
 *  HamiltonianBuilder is an abstract base class for quantum chemical methods performed with Hamiltonian parameters
 *  for which the Hamiltonian is preferably expressed as a Hermitian matrix
 *  so that the corresponding eigenvalues and -vectors can be retrieved through diagonalisation of this matrix.
 */
class HamiltonianBuilder {
public:
    // CONSTRUCTOR
    HamiltonianBuilder() = default;


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~HamiltonianBuilder() = 0;


    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @return the Hamiltonian matrix as an Eigen::MatrixXd given @param hamiltonian_parameters
     */
    virtual Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) = 0;

    /**
     *  @return the action of the Hamiltonian (@param hamiltonian_parameters and @param diagonal) on the coefficient vector @param x
     */
    virtual Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) = 0;

    /**
     *  @return the diagonal of the matrix representation of the Hamiltonian given @param hamiltonian_parameters
     */
    virtual Eigen::VectorXd calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) = 0;

    /**
     *  @return the fock space of the HamiltonianBuilder
     */
    virtual BaseFockSpace* get_fock_space() = 0;
};


}  // namespace GQCP



#endif  // GQCP_HAMILTONIANBUILDER_HPP
