#ifndef GQCG_FCI_HPP
#define GQCG_FCI_HPP


#include "HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "FockSpace/FockSpaceProduct.hpp"



namespace GQCG {

/**
 *  Full configuration interaction builds a hamiltonian matrix
 *  based on a wavefunction containing all configuration pertaining to a fixed amount of alpha and beta electrons.
 *  This means that the ONV is a combination of two ONVs, one from an alpha and beta Fock space.
 */
class FCI : public GQCG::HamiltonianBuilder {
private:
    FockSpaceProduct fock_space;  // fock space containing the alpha and beta Fock space

    struct OneElectronCoupling {
        int sign;
        size_t p;
        size_t q;
        size_t address;
    };

    // The following are rectangular arrays of dimension (dim_alpha * N_alpha * (K + 1 - N_alpha)) and similarly for beta,
    // storing one-electron excited coupling addresses (cfr. the documentation about the OneElectronCoupling struct)
    std::vector<std::vector<OneElectronCoupling>> alpha_one_electron_couplings;
    std::vector<std::vector<OneElectronCoupling>> beta_one_electron_couplings;


public:

    // CONSTRUCTORS
    /**
     *  Constructor given a @param fock_space
     */
    explicit FCI(FockSpaceProduct fock_space);


    // DESTRUCTOR
    ~FCI() = default;


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return Hamiltonian matrix as an Eigen::MatrixXd given @param hamiltonian_parameters
     */
    Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) override;

    /**
     *  @return the action of the Hamiltonian (@param hamiltonian_parameters and @param diagonal) on the coefficient vector @param x
     */
    Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) override;

    /**
     *  @return the diagonal of the matrix representation of the Hamiltonian given @param hamiltonian_parameters
     */
    Eigen::VectorXd calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) override;

    /**
     *  @return the fock space of the HamiltonianBuilder
     */
    BaseFockSpace* get_fock_space() override { return &fock_space; }
};


}  // namespace GQCG


#endif //GQCG_FCI_HPP
