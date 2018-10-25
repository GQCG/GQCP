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
#include "FockSpace/FockSpaceProduct.hpp"



namespace GQCP {


using PassToMethod = std::function<void (size_t I, size_t J, double value)>;


/**
 *  Hubbard builds a hamiltonian matrix
 *  based on a wavefunction containing all configurations pertaining to a fixed number of alpha and beta electrons.
 *  This means that a total ONV would be a combination of two ONVs, one from an alpha and one from a beta Fock space.
 *
 *  Hubbard distinguishes itself from full configuration interaction by using a simplified hamiltonian model
 *  For the two electron operators only on-site (doubly occupied in-place) interactions are considered
 *  For the one electron operators only inter-site interactions are considered
 */
class Hubbard : public GQCP::HamiltonianBuilder {
private:
    FockSpaceProduct fock_space;  // fock space containing the alpha and beta Fock space

    
    // PRIVATE METHODS
    /**
     *  Private member that evaluates either the one-electron operators for alpha or beta given the parameters.
     *  Additionally stores this evaluation in either the matvec or matrix depending on passed 
     *  @param fock_space_target refers to which spin function will be evaluated
     *  @param fock_space_fixed refers to which spin function is not evaluated
     *  @param target_is_major refers to whether or not alpha is the evaluated spin function
     *  @param hamiltonian_parameters contains the data of the operators
     *  @param method refers to which method is used (matvec or Hamliltonian matrix)
     */
    void oneOperatorModule(FockSpace& fock_space_target, FockSpace& fock_space_fixed, bool target_is_major, const HamiltonianParameters& hamiltonian_parameters, const PassToMethod& method);


public:

    // CONSTRUCTORS
    /**
     *  Constructor given a @param fock_space
     */
    explicit Hubbard(const FockSpaceProduct& fock_space);


    // DESTRUCTOR
    ~Hubbard() = default;


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return the Hamiltonian matrix as an Eigen::MatrixXd given @param hamiltonian_parameters
     */
    Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) override;

    /**
     *  @return the action of the Hamiltonian (@param hamiltonian_parameters and @param diagonal) on the coefficient vector @param x
     */
    Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) override;

    /**
     *  @return the diagonal of the matrix representation of the Hamiltonian given @param hamiltonian_parameters
     */
    Eigen::VectorXd calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) override;

    /**
     *  @return the fock space of the HamiltonianBuilder
     */
    BaseFockSpace* get_fock_space() override { return &fock_space; }
};


//  RELEVANT (non-class) METHODS
/**
 *  Generates a upper triagonal (vector) for a hubbard lattice.
 *  @param hopping_matrix allowed interaction between sites
 *  @param t one electron hopping interaction parameter
 *  @param U two electron doubly occupied interaction parameter
 *  @return the triagonal of the matrix resulting in the recominbation of U and t with the hopping matrix.
 */
Eigen::VectorXd genrateUpperTriagonal(Eigen::MatrixXd hopping_matrix, double t, double U);


}  // namespace GQCP


#endif  // GQCP_HUBBARD_HPP
