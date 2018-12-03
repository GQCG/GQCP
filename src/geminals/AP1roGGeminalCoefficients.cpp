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
#include "geminals/AP1roGGeminalCoefficients.hpp"

#include "FockSpace/FockSpace.hpp"
#include "miscellaneous.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Default constructor setting everything to zero
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients() :
    BaseAPIGGeminalCoefficients()
{}

/**
 *  @param g        the geminal coefficients in a vector representation that is in row-major storage
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients(const Eigen::VectorXd& g, size_t N_P, size_t K) :
    BaseAPIGGeminalCoefficients(g, N_P, K)
{
    if (AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K) != g.size()) {
        throw std::invalid_argument("The specified N_P and K are not compatible with the given vector of geminal coefficients.");
    }
}


/**
 *  Constructor that sets the geminal coefficients to zero
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients(size_t N_P, size_t K) :
    AP1roGGeminalCoefficients(Eigen::VectorXd::Zero(AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K)), N_P, K)
{}


/**
 *  @param ham_par      the Hamiltonian parameters
 *  @param N_P          the number of orbitals
 *
 *  @return the AP1roG geminal coefficients in the weak interaction limit
 */
AP1roGGeminalCoefficients AP1roGGeminalCoefficients::WeakInteractionLimit(const HamiltonianParameters& ham_par, size_t N_P) {

    auto K = ham_par.get_K();
    auto number_of_geminal_coefficients = AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K);
    auto h = ham_par.get_h();  // core Hamiltonian integrals
    auto g = ham_par.get_g();  // two-electron integrals

    // Provide the weak interaction limit values for the geminal coefficients
    Eigen::VectorXd g_vector = Eigen::VectorXd::Zero(number_of_geminal_coefficients);
    for (size_t mu = 0; mu < number_of_geminal_coefficients; mu++) {
        size_t i = GQCP::matrixIndexMajor(mu, K, N_P);
        size_t a = GQCP::matrixIndexMinor(mu, K, N_P);

        g_vector(mu) = - g(a,i,a,i) / (2 * (h(a,a) - h(i,i)));
    }


    return AP1roGGeminalCoefficients(g_vector, N_P, K);
}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 *
 *  @return the number of 'free' geminal coefficients
 */
size_t AP1roGGeminalCoefficients::numberOfGeminalCoefficients(size_t N_P, size_t K) {

    // Check if we can have N_P geminals in K orbitals
    if (N_P >= K) {
        throw std::invalid_argument("Can't have that many geminals in this few number of orbitals.");
    }

    return N_P * (K - N_P);
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the geminal coefficients in matrix form
 */
Eigen::MatrixXd AP1roGGeminalCoefficients::asMatrix() const {

    // Initialize the geminal coefficient matrix
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(this->N_P, this->K);

    // The AP1roG coefficients are the identity matrix in the leftmost (N_P x N_P)-block
    G.topLeftCorner(this->N_P, this->N_P) = Eigen::MatrixXd::Identity(this->N_P, this->N_P);

    // Set the right AP1roG coefficient block
    using RowMajorMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    Eigen::RowVectorXd g_row = this->g;
    RowMajorMatrixXd B = Eigen::Map<RowMajorMatrixXd, Eigen::RowMajor>(g_row.data(), this->N_P, this->K-this->N_P);
    G.topRightCorner(this->N_P, this->K-this->N_P) = B;

    return G;
}


/**
 *  @param vector_index     the vector index of the geminal coefficient
 *
 *  @return the major (geminal, non-contiguous) index i (i.e. the subscript) in the matrix of the geminal coefficients. Note that i is in [0 ... N_P[
 */
size_t AP1roGGeminalCoefficients::matrixIndexMajor(size_t vector_index) const {

    return GQCP::matrixIndexMajor(vector_index, this->K, this->N_P);
}


/**
 *  @param vector_index     the vector index of the geminal coefficient
 *
 *  @return the minor (virtual orbital, contiguous) index a (i.e. the subscript) in the matrix of the geminal coefficients. Note that a is in [N_P ... K[
 */
size_t AP1roGGeminalCoefficients::matrixIndexMinor(size_t vector_index) const {

    return GQCP::matrixIndexMinor(vector_index, this->K, this->N_P);
}


/**
 *  @param i        the major (geminal) index (changes in i are not contiguous)
 *  @param a        the minor (virtual orbital) index (changes in a are contiguous)
 *
 *  @return the vector index of the geminal coefficient G_i^a
 */
size_t AP1roGGeminalCoefficients::vectorIndex(size_t i, size_t a) const {

    if (i >= N_P) {
        throw std::invalid_argument("The major index i (subscript) must be smaller than N_P.");
    }
    if (a < N_P) {
        throw std::invalid_argument("The minor index a (superscript) must be larger than or equal to N_P.");
    }


    return GQCP::vectorIndex(i, a, this->K, this->N_P);
}


/**
 *  @return the wave function expansion corresponding to the geminal coefficients
 */
WaveFunction AP1roGGeminalCoefficients::toWaveFunction() const {

    FockSpace fock_space (this->K, this->N_P);  // the DOCI Fock space

    Eigen::MatrixXd G = this->asMatrix();  // geminal coefficients as a matrix


    Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(fock_space.get_dimension());  // coefficient vector
    ONV onv = fock_space.get_ONV(0);  // start with address 0
    for (size_t I = 0; I < fock_space.get_dimension(); I++) {

        // Construct the matrix G(m) which only has the occupied columns of G in the ONV m
        Eigen::MatrixXd Gm = Eigen::MatrixXd::Zero(this->N_P, this->N_P);

        // TODO: wait until the syntax G(Eigen::placeholders::all, occupation_indices) is released in a stable Eigen release
        for (size_t e = 0; e < this->N_P ; e++) {  // loop over all electrons
            size_t occupation_index = onv.get_occupied_index(e);

            Gm.col(e) = G.col(occupation_index);
        }


        // Calculate the permanent of Gm to obtain the coefficient
        coefficients(I) = permanent_ryser(Gm);


        if (I < fock_space.get_dimension() - 1) {  // skip the last permutation
            fock_space.setNext(onv);
        }
    }

    return WaveFunction(fock_space, coefficients);
}


}  // namespace GQCP
