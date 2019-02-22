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
#include "RHF/RHF.hpp"


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  Default constructor setting everything to zero
 */
RHF::RHF() :
    RHF(0.0, Eigen::MatrixXd::Zero(0, 0), Eigen::VectorXd::Zero(0))
{}


/**
 *  Constructor based on given converged solutions of the RHF SCF equations
 *
 *  @param electronic_energy    the converged RHF electronic energy
 *  @param C                    the coefficient matrix, i.e. the transformation matrix from the AO basis to the RHF MO basis
 *  @param orbital_energies     the RHF MO energies
 */
RHF::RHF(double electronic_energy, const Eigen::MatrixXd& C, const Eigen::VectorXd& orbital_energies) :
    electronic_energy (electronic_energy),
    C (C),
    orbital_energies (orbital_energies)
{}



/*
 *  HELPER METHODS
 */
/**
 *  @param K    the number of spatial orbitals
 *  @param N    the number of electrons
 *
 *  @return the RHF 1-RDM expressed in an orthonormal basis
 */
OneRDM calculateRHF1RDM(size_t K, size_t N) {

    if (N % 2 != 0) {
        throw std::invalid_argument("The number of given electrons cannot be odd for RHF.");
    }

    // The 1-RDM for RHF looks like (for K=5, N=6)
    //    2  0  0  0  0
    //    0  2  0  0  0
    //    0  0  2  0  0
    //    0  0  0  0  0
    //    0  0  0  0  0

    Eigen::MatrixXd D_MO = Eigen::MatrixXd::Zero(K, K);
    D_MO.topLeftCorner(N/2, N/2) = 2 * Eigen::MatrixXd::Identity(N/2, N/2);

    return OneRDM(D_MO);
}


/**
 *  @param C    the coefficient matrix, specifying the transformation to the AO basis
 *  @param N    the number of electrons
 *
 *  @return the RHF 1-RDM expressed in the AO basis
 */
Eigen::MatrixXd calculateRHFAO1RDM(const Eigen::MatrixXd& C, size_t N) {

    size_t K = C.rows();
    Eigen::MatrixXd D_MO = calculateRHF1RDM(K, N).get_matrix_representation();

    // Transform the MO 1-RDM to an AO basis
    return C * D_MO * C.adjoint();
}


/**
 *  Calculate the RHF Fock matrix F = H_core + G, in which G is a contraction of the density matrix and the two-electron integrals
 *
 *  @param D_AO     the RHF density matrix in AO basis
 *  @param ham_par  The Hamiltonian parameters in AO basis
 *
 *  @return the RHF Fock matrix expressed in the AO basis
 */
Eigen::MatrixXd calculateRHFAOFockMatrix(const Eigen::MatrixXd& D_AO, const HamiltonianParameters<double>& ham_par) {

    // To perform the contraction, we will first have to convert the Eigen::MatrixXd D_AO to an Eigen::Tensor<const double, 2> D_AO_tensor, as contractions are only implemented for Eigen::Tensors
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> D_AO_tensor (D_AO.data(), D_AO.rows(), D_AO.cols());

    // Specify the contraction pairs
    // To calculate G, we must perform two double contractions
    //      1. (mu nu|rho lambda) P(lambda rho)
    Eigen::array<Eigen::IndexPair<int>, 2> direct_contraction_pair = {Eigen::IndexPair<int>(3, 0), Eigen::IndexPair<int>(2, 1)};
    //      2. -0.5 (mu lambda|rho nu) P(lambda rho)
    Eigen::array<Eigen::IndexPair<int>, 2> exchange_contraction_pair = {Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(2, 1)};

    // Calculate both contractions (and incorporate prefactors)
//    Eigen::Tensor<double, 4> g = ham_par.get_g().get_matrix_representation();  // two-electron integrals
    auto g = ham_par.get_g();
    Eigen::Tensor<double, 2> direct_contraction = g.contract(D_AO_tensor, direct_contraction_pair);
    Eigen::Tensor<double, 2> exchange_contraction = -0.5 * g.contract(D_AO_tensor, exchange_contraction_pair);

    // The previous contractions are Eigen::Tensor<double 2> instances. In order to calculate the total G matrix, we will convert them back into Eigen::MatrixXd
    Eigen::Map<Eigen::MatrixXd> G1 (direct_contraction.data(), direct_contraction.dimension(0), direct_contraction.dimension(1));
    Eigen::Map<Eigen::MatrixXd> G2 (exchange_contraction.data(), exchange_contraction.dimension(0), exchange_contraction.dimension(1));


    // Return the final result
//    Eigen::MatrixXd H_core = ham_par.get_h().get_matrix_representation();
    auto H_core = ham_par.get_h();

    return H_core + G1 + G2;
}


/**
 *  @param D_AO         the RHF density matrix in AO basis
 *  @param H_core_AO    the core Hamiltonian parameters in AO basis
 *  @param F_AO         the Fock matrix in AO basis
 *
 *  @return the RHF electronic energy
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

    // As the double contraction of two matrices is a scalar (a tensor of rank 0), we should access the value as (0)
    return contraction(0);
}


/**
 *  @param N    the number of electrons
 *
 *  @return the RHF HOMO index
 */
size_t RHFHOMOIndex(size_t N) {

    if (N % 2 != 0) {
        throw std::invalid_argument("Can't calculate the RHF HOMO index for an odd number of electrons N.");
    }

    return N / 2 - 1;  // need to subtract 1 because computer indices start at 0
}


/**
 *  @param K    the number of spatial orbitals
 *  @param N    the number of electrons
 *
 *  @return the RHF LUMO index
 */
size_t RHFLUMOIndex(size_t K, size_t N) {

    if (N >= 2 * K) {
        throw std::invalid_argument("There is no LUMO for the given amount of electrons N and spatial orbitals K");
    }

    return RHFHOMOIndex(N) + 1;
}


}  // namespace GQCP
