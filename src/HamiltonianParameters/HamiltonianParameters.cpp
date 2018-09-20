#include "HamiltonianParameters/HamiltonianParameters.hpp"


namespace GQCG {


/**
 *  Constructor based on a given @param ao_basis_ptr, overlap @param S, one-electron operator @param h, two-electron
 *  operator @param g and a transformation matrix between the current molecular orbitals and the atomic orbitals
 *  @param C
 */
HamiltonianParameters::HamiltonianParameters(std::shared_ptr<GQCG::AOBasis> ao_basis_sptr, const GQCG::OneElectronOperator& S, const GQCG::OneElectronOperator& h, const GQCG::TwoElectronOperator& g, const Eigen::MatrixXd& C) :
    BaseHamiltonianParameters(std::move(ao_basis_sptr)),
    S (S),
    h (h),
    g (g),
    C (C)
{
    // Check if the dimensions of all matrix representations are compatible
    auto K = this->ao_basis_sptr->number_of_basis_functions;

    if ((S.dim != K) || (h.dim != K) || (g.dim != K) || (C.cols() != K) || (C.rows() != K)) {
        throw std::invalid_argument("The dimensions of the operators and coefficient matrix are incompatible.");
    }
}



}  // namespace GQCG
