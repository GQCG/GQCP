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
    HamiltonianParameters ham_par = HamiltonianParameters();  // set ham_par with accompanied memory storage
    bool is_ham_par_set = false;

    std::vector<Eigen::SparseMatrix<double>> alpha_resolved; // one-electron pq couplings
    std::vector<Eigen::SparseMatrix<double>> beta_resolved;  // two-electron rs pairs

    Eigen::SparseMatrix<double> alpha_ev; // separated evaluations (dim_alpha, dim_alpha);
    Eigen::SparseMatrix<double> beta_ev; // (dim_beta, dim_beta);


    // PRIVATE METHODS
    /**
     *  Calculates all Hamiltonian elements for operators exclusively operating for one spin function
     *  and stores these in a sparse matrix
     *
     *  @param fock_space                   Fock space for the spin function specific Hamiltonian
     *  @param k                            Modified one-electron operator
     *  @param hamiltonian_parameters       The Hamiltonian parameters in an orthonormal orbital basis
     *  @param sparse_mat                   The representation of the spin function specific Hamiltonian
     */
    void spinSeparatedModule(FockSpace& fock_space, const OneElectronOperator& k,
                             const HamiltonianParameters& hamiltonian_parameters,
                             Eigen::SparseMatrix<double>& sparse_mat) const;

    /**
     *  Calculates all one-electron couplings for the beta Fock space
     *  and attributes two-electron integrals based on the one-electron indexes of the coupling and two fixed indexes
     *
     *  @param r                        Fixed index of two-electorn integral
     *  @param s                        Fixed index of two-electron integral
     *  @param hamiltonian_parameters   The Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return                         The sparse matrix containing the calculated two-electron integrals mapped to one-electron couplings
     */
    Eigen::SparseMatrix<double> betaTwoElectronOneElectronModule(size_t r, size_t s, const HamiltonianParameters& hamiltonian_parameters) const;

    /**
     *  Calculates all one-eletron couplings for each annihilation-creation pair in the alpha Fock space
     *  and stores them in sparse matrices for each combination
     *
     *  @return vector of sparse matrices containing the one-electron couplings for the alpha Fock space
     */
    std::vector<Eigen::SparseMatrix<double>> alphaOneElectronCouplings() const;

public:

    // CONSTRUCTORS
    /**
     *  @param fock_space       the full alpha and beta product Fock space
     */
    explicit FCI(const ProductFockSpace& fock_space);


    // DESTRUCTOR
    ~FCI() = default;


    // SETTERS
    void set_hamiltonian_parameters(const HamiltonianParameters& hamiltonian_parameters);

    
    // PUBLIC METHODS
    void clearHamiltonianParameters();
    

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
