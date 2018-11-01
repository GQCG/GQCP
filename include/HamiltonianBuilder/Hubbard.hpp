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
#ifndef GQCP_HUBBARD_HPP
#define GQCP_HUBBARD_HPP


#include "HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "FockSpace/ProductFockSpace.hpp"



namespace GQCP {


// Typedef for a function that will be passed to the evaluation, it defines to what and how the evaluated results will be added.
using PassToMethod = std::function<void (size_t I, size_t J, double value)>;


/**
 *  Hubbard builds a a Hubbard Hamiltonian matrix in the FCI Fock space
 *
 *  Hubbard distinguishes itself from FCI by explicitly implementing simplified Hamiltonian parameters:
 *      - for the one electron operators only inter-site interactions are considered
 *      - for the two electron operators only on-site (doubly occupied in-place) interactions are considered
 */
class Hubbard : public GQCP::HamiltonianBuilder {
private:
    ProductFockSpace fock_space;  // fock space containing the alpha and beta Fock space

    
    // PRIVATE METHODS
    /**
     *  Evaluate the one-electron operators for alpha or beta and store the result in a matrix-vector product or a matrix, depending on the method passed
     *
     *  @param fock_space_target        the Fock space that is used as a target, i.e. that is evaluated
     *  @param fock_space_fixed         the Fock space that is not evaluated
     *  @param target_is_major          whether or not the evaluated component is the major index
     *  @param hamiltonian_parameters   the Hubbard Hamiltonian parameters
     *  @param method                   the used method: constructHamiltonian() or matrixVectorProduct()
     */
    void oneOperatorModule(FockSpace& fock_space_target, FockSpace& fock_space_fixed, bool target_is_major, const HamiltonianParameters& hamiltonian_parameters, const PassToMethod& method);


public:

    // CONSTRUCTORS
    /**
     *  @param fock_space       the full alpha and beta product Fock space
     */
    explicit Hubbard(const ProductFockSpace& fock_space);


    // DESTRUCTOR
    ~Hubbard() = default;


    // OVERRIDDEN GETTERS
    BaseFockSpace* get_fock_space() override { return &fock_space; }


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param hamiltonian_parameters       the Hubbard Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the Hubbard Hamiltonian matrix
     */
    Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) override;

    /**
     *  @param hamiltonian_parameters       the Hubbard Hamiltonian parameters in an orthonormal orbital basis
     *  @param x                            the vector upon which the Hubbard Hamiltonian acts
     *  @param diagonal                     the diagonal of the Hubbard Hamiltonian matrix
     *
     *  @return the action of the Hubbard Hamiltonian on the coefficient vector
     */
    Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) override;

    /**
     *  @param hamiltonian_parameters       the Hubbard Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the Hubbard Hamiltonian
     */
    Eigen::VectorXd calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) override;
};



/*
 *  HELPER METHODS
 */
/**
 *  Generate the upper triagonal as a vector for a Hubbard lattice
 *
 *  @param A        the adjacency matrix that represents the allowed interaction between sites
 *  @param t        the one-electron hopping interaction parameter
 *  @param U        the two-electron interaction parameter
 *
 *  @return the upper triagonal as a vector of the hopping matrix generated from the adjacency matrix and the Hubbard parameters t and U
 */
Eigen::VectorXd generateUpperTriagonal(Eigen::MatrixXd A, double t, double U);


}  // namespace GQCP


#endif  // GQCP_HUBBARD_HPP
