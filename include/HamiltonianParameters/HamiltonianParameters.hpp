#ifndef GQCG_HAMILTONIANPARAMETERS_HPP
#define GQCG_HAMILTONIANPARAMETERS_HPP


#include "Eigen/Dense"

#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"
#include "Operator/OneElectronOperator.hpp"
#include "Operator/TwoElectronOperator.hpp"
#include "JacobiRotationParameters.hpp"



namespace GQCG {


class HamiltonianParameters : public BaseHamiltonianParameters {
private:
    OneElectronOperator S;  // overlap

    OneElectronOperator h;  // one-electron interactions (i.e. the core Hamiltonian)
    TwoElectronOperator g;  // two-electron interactions

    Eigen::MatrixXd C;  // total transformation matrix between the current (restricted) molecular orbitals and the atomic orbitals


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param ao_basis, overlap @param S, one-electron operator @param h, two-electron
     *  operator @param g and a transformation matrix between the current molecular orbitals and the atomic orbitals
     *  @param C
     */
    HamiltonianParameters(std::shared_ptr<GQCG::AOBasis> ao_basis, const GQCG::OneElectronOperator& S, const GQCG::OneElectronOperator& h, const GQCG::TwoElectronOperator& g, const Eigen::MatrixXd& C);


    /**
     *  Constructor based on given Hamiltonian parameters @param ham_par and a transformation matrix @param C.
     *
     *  If the initial Hamiltonian parameters @param ham_par are expressed in the basis B, the constructed instance represents the Hamiltonian parameters in the transformed basis B'. The basis transformation between B and B' is given by the transformation matrix @param C.
     */
    HamiltonianParameters(const GQCG::HamiltonianParameters& ham_par, const Eigen::MatrixXd& C);


    // DESTRUCTORS
    ~HamiltonianParameters() override =default;

    
    // GETTERS
    GQCG::OneElectronOperator get_h() const { return this->h; }
    GQCG::TwoElectronOperator get_g() const { return this->g; }

    
    // PUBLIC METHODS
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
    void transform(const Eigen::MatrixXd& T);

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
    void rotate(const Eigen::MatrixXd& U);

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
    void rotate(const GQCG::JacobiRotationParameters& jacobi_rotation_parameters);


    // FRIEND FUNCTIONS
    friend Eigen::MatrixXd calculateRHFAOFockMatrix(const Eigen::MatrixXd& D_AO, GQCG::HamiltonianParameters ham_par);

    // FRIEND CLASSES
    friend class RHFSCFSolver;
    friend class DIISRHFSCFSolver;
};


}  // namespace GQCG


#endif  // GQCG_HAMILTONIANPARAMETERS_HPP
