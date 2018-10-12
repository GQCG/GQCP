#ifndef GQCG_CISOLVER_HPP
#define GQCG_CISOLVER_HPP


#include <HamiltonianBuilder/HamiltonianBuilder.hpp>
#include <HamiltonianParameters/HamiltonianParameters.hpp>
#include <WaveFunction.hpp>
#include <numopt.hpp>

namespace GQCG {


/**
 *  Class which requires a HamiltonianBuilder and HamiltonianParameters
 *  so that it can find a set of eigenvalues and -vectors for the Hamiltonian (which cannot always be stored in memory)
 */
class CISolver {
private:
    std::shared_ptr<const HamiltonianBuilder> hamiltonian_builder;
    HamiltonianParameters hamiltonian_parameters;

public:
    // CONSTRUCTOR
    /**
     *  Constructor given a @param hamiltonian_builder and @param hamiltonian_parameters
     */
    CISolver(const HamiltonianBuilder& hamiltonian_builder, const HamiltonianParameters& hamiltonian_parameters);


    // PUBLIC METHODS
    /**
     *  @return WaveFunction instance from solving the Hamiltonian
     */
    WaveFunction solve(numopt::eigenproblem::BaseSolverOptions solver_options);
};


}  // namespace GQCG

#endif  // GQCG_CISOLVER_HPP
