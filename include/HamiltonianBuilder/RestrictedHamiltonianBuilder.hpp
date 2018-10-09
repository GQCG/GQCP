#ifndef GQCG_HAMILTONIANBUILDER_HPP
#define GQCG_HAMILTONIANBUILDER_HPP


#include "HamiltonianParameters/HamiltonianParameters.hpp"

#include <memory>
#include <utility>



namespace GQCG {


/**
 *  RestrictedHamiltonianBuilder is an abstract base class for quantum chemical methods performed with restricted hamiltonian parameters
 *  for which the Hamiltonian is preferably expressed as a Hermitian matrix
 *  so that the corresponding eigen-values and -vectors can be retrieved through diagonalisation of this matrix.
 */
class RestrictedHamiltonianBuilder {
protected:
    HamiltonianParameters hamiltonian_parameters;  // the hamiltonian parameters passed for the calculations

    Eigen::VectorXd diagonal;  // the diagonal of the hamiltonian matrix


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor given a @param HamiltonianParameters
     */
    explicit RestrictedHamiltonianBuilder(HamiltonianParameters hamiltonian_parameters);


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
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    virtual Eigen::VectorXd calculateDiagonal() = 0;
};


}  // namespace GQCG



#endif  // GQCG_HAMILTONIANBUILDER_HPP
