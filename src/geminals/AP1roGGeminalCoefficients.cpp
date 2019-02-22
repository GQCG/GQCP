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
#include "geminals/AP1roGGeminalCoefficients.hpp"

#include "geminals/APIGGeminalCoefficients.hpp"
#include "FockSpace/FockSpace.hpp"
#include "utilities/miscellaneous.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Default constructor setting everything to zero
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients() :
    AP1roGVariables()
{}

/**
 *  @param g        the geminal coefficients in a vector representation that is in row-major storage
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients(const Eigen::VectorXd& g, size_t N_P, size_t K) :
    AP1roGVariables(g, N_P, K)
{}


/**
 *  Constructor that sets the geminal coefficients to zero
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients(size_t N_P, size_t K) :
    AP1roGVariables(N_P, K)
{}



/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  @param ham_par      the Hamiltonian parameters
 *  @param N_P          the number of orbitals
 *
 *  @return the AP1roG geminal coefficients in the weak interaction limit
 */
AP1roGGeminalCoefficients AP1roGGeminalCoefficients::WeakInteractionLimit(const HamiltonianParameters<double>& ham_par, size_t N_P) {

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
 *  DESTRUCTOR
 */
AP1roGGeminalCoefficients::~AP1roGGeminalCoefficients() {}



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

    return AP1roGVariables::numberOfVariables(N_P, K);
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
    Eigen::RowVectorXd x_row = this->x;
    RowMajorMatrixXd B = Eigen::Map<RowMajorMatrixXd, Eigen::RowMajor>(x_row.data(), this->N_P, this->K-this->N_P);
    G.topRightCorner(this->N_P, this->K-this->N_P) = B;

    return G;
}


/**
 *  @param onv      the ONV that is being projected on
 *
 *  @return the overlap of the AP1roG wave function with the given on, i.e. the projection of the APIG wave function onto that ONV
 */
double AP1roGGeminalCoefficients::overlap(const ONV& onv) const {

    // For an AP1roG, we use a simplification for singly and doubly pair-excited ONVs


    FockSpace fock_space (this->K, this->N_P);  // the DOCI Fock space
    ONV reference = fock_space.makeONV(0);

    if (onv.countNumberOfDifferences(reference) == 0) {  // no excitations
        return 1.0;
    }

    else if (onv.countNumberOfDifferences(reference) == 2) {  // one pair excitation

        size_t i = reference.findDifferentOccupations(onv)[0];
        size_t a = onv.findDifferentOccupations(reference)[0];

        return this->operator()(i, a);
    }

    else if (onv.countNumberOfDifferences(reference) == 4) {  // two pair excitations

        auto different_occupied = reference.findDifferentOccupations(onv);
        auto different_virtual = onv.findDifferentOccupations(reference);

        size_t i = different_occupied[0];
        size_t j = different_occupied[1];
        size_t a = different_virtual[0];
        size_t b = different_virtual[1];

        return this->operator()(i, a) * this->operator()(j, b) + this->operator()(j, a) * this->operator()(i, b);
    }

    else {

        APIGGeminalCoefficients APIG (this->asMatrix());
        return APIG.overlap(onv);
    }
}


}  // namespace GQCP
