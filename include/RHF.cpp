//
//  RHF.cpp
//  gqcg
//
//  Created by Laurent Lemmens on 02/10/2018.
//  Copyright © 2018 Ghent Quantum Chemistry Group. All rights reserved.
//

#include "RHF.hpp"



namespace GQCG {

/*
 * CONSTRUCTORS
 */

RHF::RHF() :
    RHF(0.0, Eigen::MatrixXd::Zero(0, 0), Eigen::VectorXd::Zero(0))
{}


RHF::RHF(double electronic_energy, const Eigen::MatrixXd& C, const Eigen::VectorXd& orbital_energies) :
    electronic_energy (electronic_energy),
    C (C),
    orbital_energies (orbital_energies)
{}



/*
 *  METHODS
 */

/**
 *  @return the RHF 1-RDM expressed in the AO basis, given the @param coefficient matrix C and the number of electrons @param N
 */
Eigen::MatrixXd calculateRHFAO1RDM(const Eigen::MatrixXd& C, size_t N) {

    size_t K = C.rows();

    // Construct the RHF density matrix in MO basis
    Eigen::MatrixXd D_MO = Eigen::MatrixXd::Zero(K, K);
    D_MO.topLeftCorner(N/2, N/2) = 2 * Eigen::MatrixXd::Identity(N/2,N/2);

    // D_AO = C * D_MO * C^dagger
    return C * D_MO * C.adjoint();
}


/**
 *  @return the RHF Fock matrix in the AO basis, given the @param D_AO density matrix in AO basis and @param ham_par_ptr Hamiltonian parameters
 */
Eigen::MatrixXd calculateRHFAOFockMatrix(const Eigen::MatrixXd& D_AO, GQCG::HamiltonianParameters ham_par) {

    // The RHF Fock matrix is calculated as F = H + G, in which G is a contraction of the density matrix and the two-electron integrals


    // To perform the contraction, we will first have to convert the Eigen::MatrixXd D_AO to an Eigen::Tensor<const double, 2> D_AO_tensor, as contractions are only implemented for Eigen::Tensors
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> D_AO_tensor (D_AO.data(), D_AO.rows(), D_AO.cols());

    // Specify the contraction pairs
    // To calculate G, we must perform two double contractions
    //      1. (mu nu|rho lambda) P(lambda rho)
    Eigen::array<Eigen::IndexPair<int>, 2> direct_contraction_pair = {Eigen::IndexPair<int>(3, 0), Eigen::IndexPair<int>(2, 1)};
    //      2. -0.5 (mu lambda|rho nu) P(lambda rho)
    Eigen::array<Eigen::IndexPair<int>, 2> exchange_contraction_pair = {Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(2, 1)};

    // Calculate both contractions (and incorporate prefactors)
    Eigen::Tensor<double, 2> direct_contraction = ham_par.g.get_matrix_representation().contract(D_AO_tensor, direct_contraction_pair);
    Eigen::Tensor<double, 2> exchange_contraction = -0.5 * ham_par.g.get_matrix_representation().contract(D_AO_tensor, exchange_contraction_pair);

    // The previous contractions are Eigen::Tensor<double 2> instances. In order to calculate the G matrix, we will convert them back into Eigen::MatrixXd instances.
    Eigen::Map<Eigen::MatrixXd> G1 (direct_contraction.data(), direct_contraction.dimension(0), direct_contraction.dimension(1));
    Eigen::Map<Eigen::MatrixXd> G2 (exchange_contraction.data(), exchange_contraction.dimension(0), exchange_contraction.dimension(1));


    // Return the final result
    return ham_par.h.get_matrix_representation() + G1 + G2;
}


/**
 *  Calculate the RHF electronic energy based on the RHF AO density matrix @param: D_AO, the core Hamiltonian @param: H_cor_AOe and the Fock matrix @param: F_AO
 */
double calculateRHFElectronicEnergy(const Eigen::MatrixXd& D_AO, const Eigen::MatrixXd& H_core_AO, const Eigen::MatrixXd& F_AO) {

    // First, calculate the sum of H_core and F (this saves a contraction)
    Eigen::MatrixXd Z = H_core_AO + F_AO;

    // Convert the matrices Z and P to an Eigen::Tensor<double, 2> P_tensor, as contractions are only implemented for Eigen::Tensors
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> D_AO_tensor (D_AO.data(), D_AO.rows(), D_AO.cols());
    Eigen::TensorMap<Eigen::Tensor<double, 2>> Z_tensor (Z.data(), D_AO.rows(), D_AO.cols());

    // Specify the contraction pair
    // To calculate the electronic energy, we must perform a double contraction
    //      0.5 P(nu mu) Z(mu nu)
    Eigen::array<Eigen::IndexPair<int>, 2> contraction_pair = {Eigen::IndexPair<int>(0, 1), Eigen::IndexPair<int>(1, 0)};

    // Calculate the double contraction (with prefactor 0.5)
    Eigen::Tensor<double, 0> contraction = 0.5 * D_AO_tensor.contract(Z_tensor, contraction_pair);

    // As the double contraction of two matrices is a scalar (a tensor of rank 0), we can access the value as (0)
    return contraction(0);
}


}  // namespace GQCG
