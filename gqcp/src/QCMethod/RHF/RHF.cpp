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
#include "QCMethod/RHF/RHF.hpp"


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  Default constructor setting everything to zero
 */
RHF::RHF() :
    RHF(0.0, TransformationMatrix<double>::Zero(0, 0), VectorX<double>::Zero(0))
{}


/**
 *  Constructor based on given converged solutions of the RHF SCF equations
 *
 *  @param electronic_energy    the converged RHF electronic energy
 *  @param C                    the coefficient matrix, i.e. the transformation matrix from the AO basis to the RHF MO basis
 *  @param orbital_energies     the RHF MO energies
 */
RHF::RHF(double electronic_energy, const TransformationMatrix<double>& C, const VectorX<double>& orbital_energies) :
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
OneRDM<double> calculateRHF1RDM(size_t K, size_t N) {

    if (N % 2 != 0) {
        throw std::invalid_argument("calculateRHF1RDM(size_t, size_t): The number of given electrons cannot be odd for RHF.");
    }

    // The 1-RDM for RHF looks like (for K=5, N=6)
    //    2  0  0  0  0
    //    0  2  0  0  0
    //    0  0  2  0  0
    //    0  0  0  0  0
    //    0  0  0  0  0

    OneRDM<double> D_MO = OneRDM<double>::Zero(K, K);
    D_MO.topLeftCorner(N/2, N/2) = 2 * SquareMatrix<double>::Identity(N/2, N/2);

    return D_MO;
}


/**
 *  @param C    the coefficient matrix, specifying the transformation to the AO basis
 *  @param N    the number of electrons
 *
 *  @return the RHF 1-RDM expressed in the AO basis
 */
OneRDM<double> calculateRHFAO1RDM(const TransformationMatrix<double>& C, size_t N) {

    size_t K = C.rows();
    auto D_MO = calculateRHF1RDM(K, N);

    // Transform the MO 1-RDM to an AO basis
    return C.conjugate() * D_MO * C.transpose();
}


/**
 *  Calculate the RHF Fock matrix F = H_core + G, in which G is a contraction of the density matrix and the two-electron integrals
 *
 *  @param D_AO             the RHF density matrix in AO basis
 *  @param sq_hamiltonian   the Hamiltonian in an AO basis
 *
 *  @return the RHF Fock matrix expressed in the AO basis
 */
ScalarSQOneElectronOperator<double> calculateRHFAOFockMatrix(const OneRDM<double>& D_AO, const SQHamiltonian<double>& sq_hamiltonian) {

    // To perform the contraction, we will first have to convert the MatrixX<double> D_AO to an Eigen::Tensor<const double, 2> D_AO_tensor, as contractions are only implemented for Eigen::Tensors
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> D_AO_tensor (D_AO.data(), D_AO.rows(), D_AO.cols());

    // Specify the contraction pairs
    // To calculate G, we must perform two double contractions
    //      1. (mu nu|rho lambda) P(lambda rho)
    Eigen::array<Eigen::IndexPair<int>, 2> direct_contraction_pair = {Eigen::IndexPair<int>(3, 0), Eigen::IndexPair<int>(2, 1)};
    //      2. -0.5 (mu lambda|rho nu) P(lambda rho)
    Eigen::array<Eigen::IndexPair<int>, 2> exchange_contraction_pair = {Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(2, 1)};

    // Calculate both contractions (and incorporate prefactors)
    const auto& g = sq_hamiltonian.twoElectron().parameters();
    Tensor<double, 2> direct_contraction = g.contract(D_AO_tensor, direct_contraction_pair);
    Tensor<double, 2> exchange_contraction = -0.5 * g.contract(D_AO_tensor, exchange_contraction_pair);

    // The previous contractions are Eigen::Tensor<double 2> instances. In order to calculate the total G matrix, we will convert them back into MatrixX<double>
    Eigen::Map<Eigen::MatrixXd> G1 (direct_contraction.data(), direct_contraction.dimension(0), direct_contraction.dimension(1));
    Eigen::Map<Eigen::MatrixXd> G2 (exchange_contraction.data(), exchange_contraction.dimension(0), exchange_contraction.dimension(1));


    return ScalarSQOneElectronOperator<double>({sq_hamiltonian.core().parameters() + G1 + G2});
}


/**
 *  @param D_AO         the RHF density matrix in AO basis
 *  @param H_core_AO    the core Hamiltonian in an AO basis
 *  @param F_AO         the Fock matrix in AO basis
 *
 *  @return the RHF electronic energy
 */
double calculateRHFElectronicEnergy(const OneRDM<double>& D_AO, const ScalarSQOneElectronOperator<double>& H_core_AO, const ScalarSQOneElectronOperator<double>& F_AO) {

    // First, calculate the sum of H_core and F (this saves a contraction)
    ScalarSQOneElectronOperator<double> Z = H_core_AO + F_AO;

    // Convert the matrices Z and P to an Eigen::Tensor<double, 2> P_tensor, as contractions are only implemented for Eigen::Tensors
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> D_AO_tensor (D_AO.data(), D_AO.rows(), D_AO.cols());
    Eigen::TensorMap<Eigen::Tensor<double, 2>> Z_tensor (Z.parameters().data(), D_AO.rows(), D_AO.cols());

    // Specify the contraction pair
    // To calculate the electronic energy, we must perform a double contraction
    //      0.5 P(nu mu) Z(mu nu)
    Eigen::array<Eigen::IndexPair<int>, 2> contraction_pair = {Eigen::IndexPair<int>(0, 1), Eigen::IndexPair<int>(1, 0)};

    // Calculate the double contraction (with prefactor 0.5)
    Tensor<double, 0> contraction = 0.5 * D_AO_tensor.contract(Z_tensor, contraction_pair);

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
        throw std::invalid_argument("RHFHOMOIndex(size_t): Can't calculate the RHF HOMO index for an odd number of electrons N.");
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
        throw std::invalid_argument("RHFLUMOIndex(size_t, size_t): There is no LUMO for the given amount of electrons N and spatial orbitals K");
    }

    return RHFHOMOIndex(N) + 1;
}



/**
 *  Specialize the orbital Hessian for RHF
 * 
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 *  @param N_P                  the number of electron pairs
 * 
 *  @return the RHF orbital Hessian as a BlockRankFourTensor, i.e. an object with a suitable operator() implemented
 */
BlockRankFourTensor<double> calculateRHFOrbitalHessianTensor(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P) {

    const auto K = sq_hamiltonian.dimension();

    BlockRankFourTensor<double> hessian (N_P,K, 0,N_P, N_P,K, 0,N_P);  // zero-initialize an object suitable for the representation of a virtual-occupied,virtual-occupied object (ai,bj)

    // Loop over all indices (ai,bj) to construct the orbital hessian
    for (size_t a = N_P; a < K; a++) {
        for (size_t i = 0; i < N_P; i++) {
            
            for (size_t b = N_P; b < K; b++) {
                for (size_t j = 0; j < N_P; j++) {
                    hessian(a,i,b,j) = calculateRHFOrbitalHessianElement(sq_hamiltonian, N_P, a, i, b, j);
                }
            }
        }
    }

    return hessian;
}


/**
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 *  @param N_P                  the number of electron pairs
 *  @param a                    the first virtual orbital index
 *  @param i                    the first occupied orbital index
 *  @param b                    the second virtual orbital index
 *  @param j                    the second occupied orbital index
 * 
 *  @return an element of the RHF orbital Hessian
 */
double calculateRHFOrbitalHessianElement(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P, const size_t a, const size_t i, const size_t b, const size_t j) {

    const auto& g = sq_hamiltonian.twoElectron().parameters();
    double value {0.0};


    // Inactive Fock matrix part
    const auto F = sq_hamiltonian.calculateInactiveFockian(N_P).parameters();
    if (i == j) {
        value += F(a,b);
    }

    if (a == b) {
        value -= F(i,j);
    }


    // Two-electron part
    value += 4 * g(a,i,b,j) - g(a,b,i,j) - g(a,j,b,i);

    return 4 * value;
}


}  // namespace GQCP
