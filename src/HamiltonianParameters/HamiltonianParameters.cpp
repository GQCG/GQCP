#include "HamiltonianParameters/HamiltonianParameters.hpp"

#include "miscellaneous.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param ao_basis, overlap @param S, one-electron operator @param h, two-electron
 *  operator @param g and a transformation matrix between the current molecular orbitals and the atomic orbitals
 *  @param C
 */
HamiltonianParameters::HamiltonianParameters(std::shared_ptr<GQCP::AOBasis> ao_basis, const GQCP::OneElectronOperator& S, const GQCP::OneElectronOperator& h, const GQCP::TwoElectronOperator& g, const Eigen::MatrixXd& C) :
    BaseHamiltonianParameters(std::move(ao_basis)),
    K (S.dim),
    S (S),
    h (h),
    g (g),
    C (C)
{
    // Check if the dimensions of all matrix representations are compatible
    auto error = std::invalid_argument("The dimensions of the operators and coefficient matrix are incompatible.");

    if (this->ao_basis) {  // ao_basis is not nullptr
        if (this->K != this->ao_basis->number_of_basis_functions) {
            throw error;
        }
    }

    if ((h.dim != this->K) || (g.dim != this->K) || (C.cols() != this->K) || (C.rows() != this->K)) {
        throw error;
    }
}


/**
 *  Constructor based on given Hamiltonian parameters @param ham_par and a transformation matrix @param C.
 *
 *  If the initial Hamiltonian parameters @param ham_par are expressed in the basis B, the constructed instance represents the Hamiltonian parameters in the transformed basis B'. The basis transformation between B and B' is given by the transformation matrix @param C.
 */
HamiltonianParameters::HamiltonianParameters(const GQCP::HamiltonianParameters& ham_par, const Eigen::MatrixXd& C) :
    BaseHamiltonianParameters(ham_par.ao_basis),
    K (ham_par.S.dim),
    S (ham_par.S),
    h (ham_par.h),
    g (ham_par.g),
    C (ham_par.C)
{
    // We have now initialized the new Hamiltonian parameters to be a copy of the given Hamiltonian parameters, so now we will transform
    this->transform(C);
}


/*
 *  PUBLIC METHODS
 */

/**
 *  Given a transformation matrix @param T that links the new molecular orbital basis to the old molecular orbital basis,
 *  in the sense that
 *       b' = b T ,
 *  in which the molecular orbitals are collected as elements of a row vector b, transform
 *      - the one-electron interaction operator (i.e. the core Hamiltonian)
 *      - the two-electron interaction operator
 *
 *  Furthermore
 *      - @member S now gives the overlap matrix in the new molecular orbital basis
 *      - @member C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
 */
void HamiltonianParameters::transform(const Eigen::MatrixXd& T) {

    this->S.transform(T);

    this->h.transform(T);
    this->g.transform(T);

    this->C = this->C * T;  // use the correct transformation formula for subsequent transformations
}


/**
 *  Given a unitary rotation matrix @param U that links the new molecular orbital basis to the old molecular orbital basis,
 *  in the sense that
 *       b' = b U ,
 *  in which the molecular orbitals are collected as elements of a row vector b, transform
 *      - the one-electron interaction operator (i.e. the core Hamiltonian)
 *      - the two-electron interaction operator
 *
 *  Furthermore, @member C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
 */
void HamiltonianParameters::rotate(const Eigen::MatrixXd& U) {

    // A rotation leaves the overlap matrix invariant, so we don't have to transform it

    this->h.rotate(U);
    this->g.rotate(U);

    this->C = this->C * U;
}


/**
 *  Given @param jacobi_rotation_parameters that represent a unitary rotation matrix @param U (using a (cos, sin, -sin, cos) definition for the Jacobi rotation matrix) that links the new molecular orbital basis to the old molecular orbital basis,
 *  in the sense that
 *       b' = b U ,
 *  in which the molecular orbitals are collected as elements of a row vector b, transform
 *      - the one-electron interaction operator (i.e. the core Hamiltonian)
 *      - the two-electron interaction operator
 *
 *  Furthermore @member C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
 */
void HamiltonianParameters::rotate(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters) {

    // A rotation leaves the overlap matrix invariant, so we don't have to transform it

    this->h.rotate(jacobi_rotation_parameters);
    this->g.rotate(jacobi_rotation_parameters);


    // Create a Jacobi rotation matrix to transform the coefficient matrix with
    size_t K = this->h.dim;  // number of spatial orbitals
    auto J = GQCP::jacobiRotationMatrix(jacobi_rotation_parameters, K);
    this->C = this->C * J;
}


/**
 *  Given @param one_rdm and @param two_rdm
 *  @return the energy as a result of the contraction of the 1- and 2-RDMs with the one- and two-electron integrals
 */
double HamiltonianParameters::calculateEnergy(OneRDM one_rdm, TwoRDM two_rdm){

    double energy_by_contraction = (this->h.get_matrix_representation() * one_rdm.get_matrix_representation()).trace();

    Eigen::Tensor<double, 4> d = two_rdm.get_matrix_representation();
    Eigen::Tensor<double, 4> g = this->g.get_matrix_representation();

    // Specify the contractions for the relevant contraction of the two-electron integrals and the 2-RDM
    //      0.5 g(p q r s) d(p q r s)
    Eigen::array<Eigen::IndexPair<int>, 4> contractions = {Eigen::IndexPair<int>(0,0), Eigen::IndexPair<int>(1,1), Eigen::IndexPair<int>(2,2), Eigen::IndexPair<int>(3,3)};
    //      Perform the contraction
    Eigen::Tensor<double, 0> contraction = 0.5 * g.contract(d, contractions);

    // As the contraction is a scalar (a tensor of rank 0), we should access by (0).
    energy_by_contraction += contraction(0);

    return energy_by_contraction;
}

}  // namespace GQCP
