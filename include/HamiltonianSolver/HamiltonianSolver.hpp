#ifndef GQCG_HAMILTONIANSOLVER_HPP
#define GQCG_HAMILTONIANSOLVER_HPP


namespace GQCG {

/**
 *  Class which requires a HamiltonianBuilder and HamiltonianParameters
 *  so that it can be solved to find a set of eigenvalues and -vectors for the Hamiltonian (which cannot always be stored in memory)
 */
class HamiltonianSolver {
private:
    HamiltonianBuilder hamiltonian_builder;
    HamiltonianParameters hamiltonian_parameters;

public:
    // CONSTRUCTOR
    /**
     *  Constructor given a @param hamiltonian_builder and @param hamiltonian_parameters
     */
    HamiltonianSolver(HamiltonianBuilder& hamiltonian_builder, HamiltonianParameters& hamiltonian_parameters);


    // PUBLIC METHODS
    /**
     *  @return WaveFunction instance from solving the Hamiltonian
     */
    WaveFunction solve(numopt::eigenproblem::BaseSolverOptions solver_options);
};


}  // namespace GQCG

#endif  // GQCG_HAMILTONIANSOLVER_HPP
