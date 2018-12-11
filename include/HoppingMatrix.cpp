#include "HoppingMatrix.hpp"

#include "miscellaneous.hpp"
#include <iostream>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param H        the Hubbard hopping matrix
 */
HoppingMatrix::HoppingMatrix(const Eigen::MatrixXd& H) :
    K (H.cols()),
    H (H)
{
    if (H.cols() != H.rows()) {
        throw std::invalid_argument("The given hopping matrix must be square.");
    }

    if (!H.transpose().isApprox(H)) {
        throw std::invalid_argument("The given hopping matrix must be symmetric.");
    }
}



/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  @param upper_triangle       the upper triangle (in column-major ordering) that specifies the Hubbard hopping matrix
 *
 *  @return the hopping matrix that corresponds to the given upper triangle
 */
HoppingMatrix HoppingMatrix::FromUpperTriangle(const Eigen::VectorXd& upper_triangle) {

    // Check if the given upper triangle is valid
    size_t x = upper_triangle.size();
    size_t K = (static_cast<size_t>(std::sqrt(1 + 8*x) - 1))/2;  // number of rows and columns

    if (K * (K+1) != 2*x) {
        throw std::invalid_argument("The given upper triangle is not a valid upper triangle");
    }


    return HoppingMatrix(fromUpperTriangle(upper_triangle));
}


/**
 *  @param K        the number of lattice sites
 *
 *  @return a random hopping matrix with elements distributed uniformly in [-1.0, 1.0]
 */
static HoppingMatrix Random(size_t K) {

    Eigen::VectorXd v = Eigen::VectorXd::Random(K*(K+1)/2);  // random free variables

    return HoppingMatrix::FromUpperTriangle(v);
}



}  // namespace GQCP
