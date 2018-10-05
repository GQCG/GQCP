#ifndef GQCG_DOCI_HPP
#define GQCG_DOCI_HPP


#include <memory>
#include "RestrictedHamiltonianBuilder.hpp"
#include "FockSpace/FockSpace.hpp"



namespace GQCG {


class DOCI : public GQCG::RestrictedHamiltonianBuilder {
private:
    FockSpace* fock_space;
    size_t dim;


public:
    // CONSTRUCTORS
    /**
     *  constructor given a @param HamiltonianParameters and FockSpace
     */
    explicit DOCI(HamiltonianParameters &hamiltonian_parameters, FockSpace &fock_space);

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
    Eigen::VectorXd matrixVectorProduct(const Eigen::VectorXd& x)override;

    /**
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    Eigen::VectorXd calculateDiagonal() override;
};


}  // namespace GQCG


#endif //GQCG_DOCI_HPP
