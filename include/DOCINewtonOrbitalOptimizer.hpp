#ifndef DOCINewtonOrbitalOptimizer_hpp
#define DOCINewtonOrbitalOptimizer_hpp

#include "FockSpace/FockSpace.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "CISolver/CISolver.hpp"

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
    GQCG::DOCI doci;

    GQCG::HamiltonianParameters ham_par;
    numopt::eigenproblem::BaseSolverOptions& solver_options;
    GQCG::CISolver doci_solver;

    const double oo_convergence_threshold;
    const size_t maximum_number_of_oo_iterations;
    bool is_converged = false;


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param fock_space, Hamiltonian parameters @param ham_par, @param solver_options for solving the DOCI eigenvalue problem, a @param oo_convergence_threshold and a @param maximum_number_of_oo_iterations
     */
    DOCINewtonOrbitalOptimizer(const GQCG::FockSpace& fock_space, const GQCG::HamiltonianParameters& ham_par, numopt::eigenproblem::BaseSolverOptions& solver_options, double oo_convergence_threshold=1.0e-08, size_t maximum_number_of_oo_iterations=128);

    // GETTERS
    std::vector<numopt::eigenproblem::Eigenpair> get_eigenpairs() const;
    numopt::eigenproblem::Eigenpair get_eigenpair(size_t index = 0) const;

    // PUBLIC METHODS
    /**
     *  Perform the orbital optimization
     */
    void solve();

    /**
     *  @return a WaveFunction instance after performing the orbital optimization for a given eigenvector at @param index
     */
    GQCG::WaveFunction get_wavefunction(size_t index = 0);
};


}  // namespace GQCG


#endif /* DOCINewtonOrbitalOptimizer_hpp */
