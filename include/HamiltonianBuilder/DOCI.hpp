#ifndef GQCG_DOCI_HPP
#define GQCG_DOCI_HPP


#include "RestrictedHamiltonianBuilder.hpp"
#include "FockSpace/FockSpace.hpp"

#include <memory>



namespace GQCG {


/**
 *  Doubly occupied configuration interaction builds a hamiltonian matrix
 *  based on a wavefunction only containing doubly occupied configurations.
 *  This means that the combined ONV from both the alpha and beta Fock space
 *  requires the individual ONVs to be identical (beta configuration = alpha configuration).
 *  In turn this is only possible when both Fock spaces are identical.
 */
class DOCI : public GQCG::RestrictedHamiltonianBuilder {
private:
    FockSpace fock_space;  // both the alpha and beta Fock space
    size_t dim;  // dimension of this->fock_space


public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param hamiltonian_parameters and @param fock_space
     */
    explicit DOCI(HamiltonianParameters hamiltonian_parameters, FockSpace fock_space);


    // DESTRUCTOR
    ~DOCI() = default;


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return Hamiltonian matrix as an Eigen::MatrixXd
     */
    Eigen::MatrixXd constructHamiltonian() override;

    /**
     *  @return the action of the Hamiltonian of the coefficient vector @param x
     */
    Eigen::VectorXd matrixVectorProduct(const Eigen::VectorXd& x) override;

    /**
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    Eigen::VectorXd calculateDiagonal() override;
};


}  // namespace GQCG


#endif  // GQCG_DOCI_HPP
