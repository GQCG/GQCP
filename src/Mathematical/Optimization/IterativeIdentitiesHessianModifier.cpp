#include "Mathematical/Optimization/IterativeIdentitiesHessianModifier.hpp"

#include "Eigen/Cholesky"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param alpha                the increment factor used to obtain the next scaler tau
 *  @param beta                 the heuristic increment
 */
IterativeIdentitiesHessianModifier::IterativeIdentitiesHessianModifier(const double alpha, const double beta) :
    alpha (alpha),
    beta (beta)
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param hessian      the current indefinite Hessian
 * 
 *  @return a modified Hessian that is made positive (for minimizers) or negative (for maximizers) definite
 */
SquareMatrix<double> IterativeIdentitiesHessianModifier::operator()(const SquareMatrix<double>& hessian) {

    // Initialize tau_0, i.e. the initial value for tau
    const double min_diagonal_el = hessian.diagonal().minCoeff();
    if (min_diagonal_el > 0) {
        this->tau = 0.0;
    } else {
        this->tau = - min_diagonal_el + this->beta;
    }


    // Update the modified hessian with multiples of the identity matrix until it is positive/negative definite
    const size_t dim = hessian.get_dim();
    SquareMatrix<double> modified_hessian = hessian;
    Eigen::LLT<Eigen::MatrixXd> llt_factorizer {};  // use Cholesky decomposition to check for positive/negative definiteness
    while ((llt_factorizer.compute(modified_hessian), llt_factorizer.info() != Eigen::Success)) {  // comma operator gets value of last expression

        // Modify the Hessian with a multiple of the identity matrix
        modified_hessian += this->tau * SquareMatrix<double>::Identity(dim, dim);

        // Update the scaler
        this->tau = std::max(this->alpha * this->tau, this->beta);
    }

    return modified_hessian;
}


}  // namespace GQCP
