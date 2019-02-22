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
#ifndef GQCP_SELECTEDCI_HPP
#define GQCP_SELECTEDCI_HPP


#include "HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "FockSpace/SelectedFockSpace.hpp"



namespace GQCP {

/**
 *  SelectedCI builds a Hamiltonian matrix in the Selected Fock space
 */
class SelectedCI : public HamiltonianBuilder {
private:
    SelectedFockSpace fock_space;  // contains both the alpha and beta Fock space
    
    // PRIVATE METHODS
    /**
     *  Evaluate all Hamiltonian elements, putting the results in the Hamiltonian matrix or matvec through the `method` function
     *  This function is used both in `constructHamiltonian()` and `matrixVectorProduct()` to avoid duplicate code.
     *
     *  @param hamiltonian_parameters   the Hamiltonian parameters in an orthonormal basis
     *  @param method                   the method depending to how you wish to construct the Hamiltonian
     */
    void evaluateHamiltonianElements(const HamiltonianParameters<double>& hamiltonian_parameters, const PassToMethod& method) const;
public:

    // CONSTRUCTORS
    /**
     *  @param fock_space               the selected Fock space
     */
    explicit SelectedCI(const SelectedFockSpace& fock_space);


    // DESTRUCTOR
    ~SelectedCI() = default;


    // OVERRIDDEN GETTERS
    const BaseFockSpace* get_fock_space() const override { return &fock_space; }


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the SelectedCI Hamiltonian matrix
     */
    Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters<double>& hamiltonian_parameters) const override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *  @param x                            the vector upon which the SelectedCI Hamiltonian acts
     *  @param diagonal                     the diagonal of the SelectedCI Hamiltonian matrix
     *
     *  @return the action of the SelectedCI Hamiltonian on the coefficient vector
     */
    Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters<double>& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) const override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the SelectedCI Hamiltonian
     */
    Eigen::VectorXd calculateDiagonal(const HamiltonianParameters<double>& hamiltonian_parameters) const override;
};



}  // namespace GQCP


#endif  // GQCP_SELECTEDCI_HPP
