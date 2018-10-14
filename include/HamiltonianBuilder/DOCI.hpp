#ifndef GQCG_DOCI_HPP
#define GQCG_DOCI_HPP


#include "HamiltonianBuilder.hpp"
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
class DOCI : public GQCG::HamiltonianBuilder {
private:
    FockSpace fock_space;  // both the alpha and beta Fock space
    size_t dim;  // dimension of this->fock_space


public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param hamiltonian_parameters and @param fock_space
     */
    explicit DOCI(FockSpace fock_space);


    // DESTRUCTOR
    ~DOCI() = default;


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return Hamiltonian matrix as an Eigen::MatrixXd
     */
    Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) override;

    /**
     *  @return the action of the Hamiltonian of the coefficient vector @param x
     */
    Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::MatrixXd& diagonal) override;

    /**
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    Eigen::VectorXd calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) override;

    /**
     *  @return the fock space of the HamiltonianBuilder
     */
    BaseFockSpace* get_fock_space() override { return &fock_space; }
};


}  // namespace GQCG


#endif  // GQCG_DOCI_HPP
