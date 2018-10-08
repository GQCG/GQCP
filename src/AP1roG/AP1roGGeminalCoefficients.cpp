//
//  AP1roGGeminalCoefficients.cpp
//  gqcg
//
//  Created by Laurent Lemmens on 05/10/2018.
//  Copyright © 2018 Ghent Quantum Chemistry Group. All rights reserved.
//

#include "AP1roGGeminalCoefficients.hpp"
#include <iostream>

namespace GQCG {

/*
 *  CONSTRUCTORS
 */
/**
 *  Default constructor setting the geminal coefficients to zero, based on the number of orbitals @param K and number of electron pairs @param N_P
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients(size_t N_P, size_t K) :
    N_P (N_P),
    K (K),
    g (Eigen::VectorXd::Zero((K-N_P) * N_P))
{
    // Check if we can create N_P geminals in K orbitals
    if (N_P >= K) {
        throw std::invalid_argument("Can't create that many geminals in this few number of orbitals.");
    }
}


/**
 *  Constructor based on given geminal coefficients @param g, which are in vector (row-major) form
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients(const Eigen::VectorXd& g, size_t N_P, size_t K) :
    N_P (N_P),
    K (K),
    g (g)
{
    if (N_P * (K-N_P) != g.size()) {
        throw std::invalid_argument("The specified N_P and K are not compatible with the given vector of geminal coefficients.");
    }
}


/*
 *  PUBLIC METHODS
 */
/**
 *  Construct and @return the geminal coefficients in matrix form
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

}  // namespace GQCG
