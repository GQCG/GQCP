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



namespace GQCP {

/**
 *  A HamiltonianBuilder for FCI: it builds the matrix representation of the FCI Hamiltonian in the full alpha and beta product Fock space
 */
class FCI : public GQCP::HamiltonianBuilder {
private:
    ProductFockSpace fock_space;  // fock space containing the alpha and beta Fock space

    // Rectangular matrix of SpinEvaluations
    /**
     *  A small struct that is used to hold in memory the addresses of spin strings differing in one electron
     *  excitation (an annihilation on orbital p and a creation on orbital q) that are coupled through the Hamiltonian
     *
     *  During the construction of the FCI Hamiltonian, the one-electron excited coupling strings are both needed in the
     *  alpha, beta, and alpha-beta parts. When a spin string is found that couples to another spin string (with address
     *  I), the address of the coupling spin string is hold in memory, in the following way: in a
     *  std::vector<std::vector<OneElectronCoupling>> (with dimension I_alpha * N_alpha * (K + 1 - N_alpha)), at every outer index
     *  I_alpha, a std::vector of OneElectronCouplings is kept, each coupling through the Hamiltonian to that particular
     *  spin string with address I_alpha. The beta case is similar. The sign of the matrix element, i.e. <I_alpha | H | address> is also stored.
     *
     *  We can keep this many addresses in memory because the resulting dimension (cfr. dim_alpha * N_alpha * (K + 1 - N_alpha)) is
     *  significantly less than the dimension of the FCI space (cfr. I_alpha * I_beta).
     *
     *  The number of coupling spin strings for an alpha string is equal to N_alpha * (K + 1 - N_alpha), since we have to pick
     *  one out of N_alpha occupied indices to annihilate, and afterwards (after the annihilation) we have (K + 1 - N_A)
     *  choices to pick an index to create on.
     */
    struct OneElectronCoupling {
        int sign;
        size_t p;
        size_t q;
        size_t address;
    };

    // The following are rectangular arrays of dimension (dim_alpha * N_alpha * (K + 1 - N_alpha)) and similarly for beta,
    // storing one-electron excited coupling addresses (cfr. the documentation about the OneElectronCoupling struct)
    std::vector<std::vector<OneElectronCoupling>> alpha_one_electron_couplings;
    std::vector<std::vector<OneElectronCoupling>> beta_one_electron_couplings;


public:

    // CONSTRUCTORS
    /**
     *  @param fock_space       the full alpha and beta product Fock space
     */
    explicit FCI(const ProductFockSpace& fock_space);


    // DESTRUCTOR
    ~FCI() = default;


    // OVERRIDDEN GETTERS
    BaseFockSpace* get_fock_space() override { return &fock_space; }


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the FCI Hamiltonian matrix
     */
    Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *  @param x                            the vector upon which the FCI Hamiltonian acts
     *  @param diagonal                     the diagonal of the FCI Hamiltonian matrix
     *
     *  @return the action of the FCI Hamiltonian on the coefficient vector
     */
    Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    Eigen::VectorXd calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) override;
};


}  // namespace GQCP


#endif //GQCP_FCI_HPP
