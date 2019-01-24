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
#ifndef GQCP_FCI_HPP
#define GQCP_FCI_HPP


#include "HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "FockSpace/ProductFockSpace.hpp"

#include <Eigen/Sparse>


namespace GQCP {

/**
 *  A HamiltonianBuilder for FCI: it builds the matrix representation of the FCI Hamiltonian in the full alpha and beta product Fock space
 */
class FCI : public HamiltonianBuilder {
private:
    ProductFockSpace fock_space;  // fock space containing the alpha and beta Fock space
    std::vector<Eigen::SparseMatrix<double>> alpha_couplings;

    // PRIVATE METHODS
    /**
     *  Calculates the matrix Hamiltonian in the alpha or beta Fock space
     *
     *  @param fock_space                   Fock space for the spin function specific Hamiltonian
     *  @param hamiltonian_parameters       The Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return The sparse matrix containing all Hamiltonian elements for the Fock space pertaining to a single spin
     */
    Eigen::SparseMatrix<double> calculateSpinSeparatedHamiltonian(const FockSpace& fock_space, const HamiltonianParameters& hamiltonian_parameters) const;

    /**
     *  Calculates theta[rs]: all one-electron couplings for a (spin) Fock space
     *  and attributes two-electron integrals based on the one-electron coupling and two chosen fixed indexes
     *
     *  @param r                        First index of the two-electron integral
     *  @param s                        Second index of the two-electron integral
     *  @param fock_space
     *  @param hamiltonian_parameters   The Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return The sparse matrix containing the calculated two-electron integrals mapped to one-electron couplings
     */
    Eigen::SparseMatrix<double> calculateTwoElectronIntermediate(size_t r, size_t s, const HamiltonianParameters& hamiltonian_parameters, const FockSpace& fock_space) const;

    /**
     *  Calculates sigma(pq) + sigma(qp)'s: all one-electron couplings for each annihilation-creation pair in the (spin) Fock space
     *  and stores them in sparse matrices for each combination pair
     *
     *  @return vector of sparse matrices containing the one-electron couplings for the (spin) Fock space
     *  Ordered as: sigma(00), sigma(01) + sigma(10), sigma(02)+ sigma(20), ...
     */
    std::vector<Eigen::SparseMatrix<double>> calculateOneElectronCouplingsIntermediates(const FockSpace& fock_space) const;
public:

    // CONSTRUCTORS
    /**
     *  @param fock_space       the full alpha and beta product Fock space
     */
    explicit FCI(const ProductFockSpace& fock_space);


    // DESTRUCTOR
    ~FCI() = default;


    // OVERRIDDEN GETTERS
    const BaseFockSpace* get_fock_space() const override { return &fock_space; }


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the FCI Hamiltonian matrix
     */
    Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) const override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *  @param x                            the vector upon which the FCI Hamiltonian acts
     *  @param diagonal                     the diagonal of the FCI Hamiltonian matrix
     *
     *  @return the action of the FCI Hamiltonian on the coefficient vector
     */
    Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) const override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    Eigen::VectorXd calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) const override;
};


}  // namespace GQCP


#endif //GQCP_FCI_HPP
