#ifndef GQCG_CISOLVER_HPP
#define GQCG_CISOLVER_HPP


#include "HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "WaveFunction.hpp"

#include <numopt.hpp>


namespace GQCG {


/**
 *  Class which solves the CI eigenvalue problem and requires a HamiltonianBuilder and HamiltonianParameters
 *  so that it can find a set of eigenvalues and -vectors for the Hamiltonian (which cannot always be stored in memory)
 */
class CISolver {
private:
    HamiltonianBuilder* hamiltonian_builder;
    HamiltonianParameters hamiltonian_parameters;

    std::vector<numopt::eigenproblem::Eigenpair> eigenpairs;  // eigenvalues and -vectors

public:
    // CONSTRUCTOR
    /**
     *  Constructor given a @param hamiltonian_builder and @param hamiltonian_parameters
     */
    CISolver(HamiltonianBuilder& hamiltonian_builder, const HamiltonianParameters& hamiltonian_parameters);


    // GETTERS
    std::vector<numopt::eigenproblem::Eigenpair> get_eigenpairs() const { return this->eigenpairs; }
    numopt::eigenproblem::Eigenpair get_eigenpair(size_t index = 0) const { return this->eigenpairs[index]; }


    // PUBLIC METHODS
    /**
     *  solves the CI problem, setting the eigenpairs
     */
    void solve(numopt::eigenproblem::BaseSolverOptions& solver_options);

    /**
     *  @return WaveFunction instance after solving the CI problem for a given eigenvector at @param index
     */
    GQCG::WaveFunction get_wavefunction(size_t index = 0);
};


}  // namespace GQCG

#endif  // GQCG_CISOLVER_HPP
