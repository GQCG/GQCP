#include <utility>

#ifndef GQCG_HAMILTONIANBUILDER_HPP
#define GQCG_HAMILTONIANBUILDER_HPP


#include <memory>
#include "HamiltonianParameters/HamiltonianParameters.hpp"



namespace GQCG {


class RestrictedHamiltonianBuilder {
protected:
    OneElectronOperator h;
    TwoElectronOperator g;
    HamiltonianParameters hamiltonian_parameters;  // the hamiltonian parameters passed for the calculations

    Eigen::VectorXd diagonal;


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor given a @param HamiltonianParameters
     */
    explicit RestrictedHamiltonianBuilder(HamiltonianParameters hamiltonian_parameters) :
    hamiltonian_parameters(hamiltonian_parameters),
    h(hamiltonian_parameters.h),
    g(hamiltonian_parameters.g) {}


public:
    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~RestrictedHamiltonianBuilder() = 0;


    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @return Hamiltonian matrix as an Eigen::MatrixXd
     */
    virtual Eigen::MatrixXd constructHamiltonian() = 0;

    /**
     *  @return the action of the Hamiltonian of the coefficient vector @param x
     */
    virtual Eigen::VectorXd matrixVectorProduct(const Eigen::VectorXd& x) = 0;

    /**
     *  @set the diagonal of the matrix representation of the Hamiltonian
     */
    virtual Eigen::VectorXd calculateDiagonal() = 0;
};


}  // namespace GQCG



#endif //GQCG_HAMILTONIANBUILDER_HPP
