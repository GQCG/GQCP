#ifndef DOCINewtonOrbitalOptimizer_hpp
#define DOCINewtonOrbitalOptimizer_hpp

#include "FockSpace/FockSpace.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"

#include <numopt.hpp>


namespace GQCG {


/**
 *  A class that performs gradient-and-Hessian-based orbital optimization for DOCI by sequentially
 *      - solving the DOCI eigenvalue problem
 *      - solving the Newton step to find the anti-Hermitian orbital rotation parameters
 *      - rotating the underlying spatial orbital basis
 */
class DOCINewtonOrbitalOptimizer {
private:
    const GQCG::FockSpace fock_space;

    GQCG::HamiltonianParameters ham_par;
    numopt::eigenproblem::BaseSolverOptions& solver_options;
    const double oo_convergence_threshold;
    const size_t maximum_number_of_oo_iterations;
    bool is_converged = false;


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param fock_space, Hamiltonian parameters @param ham_par, @param solver_options for solving the DOCI eigenvalue problem, a @param oo_convergence_threshold and a @param maximum_number_of_oo_iterations
     */
    DOCINewtonOrbitalOptimizer(const GQCG::FockSpace& fock_space, const GQCG::HamiltonianParameters& ham_par, numopt::eigenproblem::BaseSolverOptions& solver_options, double oo_convergence_threshold=1.0e-08, size_t maximum_number_of_oo_iterations=128);


    // PUBLIC METHODS
    /**
     *  Perform the orbital optimization
     */
    void solve();

};


}  // namespace GQCG


#endif /* DOCINewtonOrbitalOptimizer_hpp */
