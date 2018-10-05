#ifndef RHF_hpp
#define RHF_hpp


#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


namespace GQCG {

/**
 *  A class that represents a converged solution to the RHF SCF equations. It has
 *      - @member electronic_energy: the converged RHF electronic energy
 *      - @member C: the coefficient matrix, i.e. the matrix that links the AO basis to the RHF MO basis
 *      - @member orbital_energies: the energies of the RHF MOs
 */
class RHF {
private:
    double electronic_energy;
    Eigen::MatrixXd C;  // transformation matrix from the AO basis to the RHF MO basis
    Eigen::VectorXd orbital_energies;  // sorted in ascending energies




public:
    // CONSTRUCTORS
    /**
     *  Default constructor setting everything to zero
     */
    RHF();

    /**
     *  Constructor based on given converged solutions of the RHF SCF equations
     */
    RHF(double electronic_energy, const Eigen::MatrixXd& C, const Eigen::VectorXd& orbital_energies);


    // GETTERS
    double get_electronic_energy() const { return this->electronic_energy; }
    Eigen::MatrixXd get_C() const { return this->C; }
    Eigen::VectorXd get_orbital_energies() const { return this->orbital_energies; }

    // FRIEND CLASSES
    friend class RHFSCFSolver;
};


/*
 *  HELPER METHODS
 */
/**
 *  @return the RHF 1-RDM expressed in the AO basis, given the @param coefficient matrix C and the number of electrons @param N
 */
Eigen::MatrixXd calculateRHFAO1RDM(const Eigen::MatrixXd& C, size_t N);

/**
 *  @return the RHF Fock matrix in the AO basis, given the @param D_AO density matrix in AO basis and @param ham_par_ptr Hamiltonian parameters
 *
 *  The RHF Fock matrix is calculated as F = H + G, in which G is a contraction of the density matrix and the two-electron integrals
 */
Eigen::MatrixXd calculateRHFAOFockMatrix(const Eigen::MatrixXd& D_AO, GQCG::HamiltonianParameters ham_par);

/**
 *  @return the RHF electronic energy based on the RHF AO density matrix @param: D_AO, the core Hamiltonian @param: H_core_AO and the Fock matrix @param: F_AO
 */
double calculateRHFElectronicEnergy(const Eigen::MatrixXd& D_AO, const Eigen::MatrixXd& H_core_AO, const Eigen::MatrixXd& F_AO);

}  // namespace GQCG


#endif /* RHF_hpp */
