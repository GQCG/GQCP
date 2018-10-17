#ifndef GQCG_FCI_HPP
#define GQCG_FCI_HPP


#include "HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "FockSpace/FockSpaceProduct.hpp"



namespace GQCG {

/**
 *  Full configuration interaction (FCI) builds a hamiltonian matrix
 *  based on a wavefunction containing all configurations pertaining to a fixed number of alpha and beta electrons.
 *  This means that a total ONV would be a combination of two ONVs, one from an alpha and one from a beta Fock space.
 */
class FCI : public GQCG::HamiltonianBuilder {
private:
    FockSpaceProduct fock_space;  // fock space containing the alpha and beta Fock space

    // Rectangular matrix of SpinEvaluations
    /**
     *  A small struct that is used to hold in memory the @param address of spin strings differing in one electron
     *  excitation (an annihilation on orbital @param p and a creation on orbital @param q) that are coupled through the
     *  Hamiltonian.
     *
     *  During the construction of the FCI Hamiltonian, the one-electron excited coupling strings are both needed in the
     *  alpha, beta, and alpha-beta parts. When a spin string is found that couples to another spin string (with address
     *  I), the address of the coupling spin string is hold in memory, in the following way: in a
     *  std::vector<std::vector<OneElectronCoupling>> (with dimension I_alpha * N_alpha * (K + 1 - N_alpha)), at every outer index
     *  I_alpha, a std::vector of OneElectronCouplings is kept, each coupling through the Hamiltonian to that particular
     *  spin string with address I_alpha. Of course, the beta case is similar.
     *
     *  The @param sign of the matrix element, i.e. <I_alpha | H | address> is also stored as a parameter.
     *
     *
     *  We can keep this many addresses in memory because the resulting dimension (cfr. dim_alpha * N_alpha * (K + 1 - N_alpha)) is
     *  significantly less than the dimension of the FCI space (cfr. I_alpha * I_beta).
     *
     *  The number of coupling spin strings for an alpha string is equal to N_alpha * (K + 1 - N_alpha), since we have to pick
     *  one out of N_alpha occupied indices to annihilate, and afterwards (after the annihilation) we have (K + 1 - N_A)
     *  choices to pick an index to create on.
     */
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
    explicit FCI(const FockSpaceProduct& fock_space);


    // DESTRUCTOR
    ~FCI() = default;


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return the Hamiltonian matrix as an Eigen::MatrixXd given @param hamiltonian_parameters
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
