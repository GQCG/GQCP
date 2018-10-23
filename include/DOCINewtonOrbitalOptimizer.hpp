#ifndef DOCINewtonOrbitalOptimizer_hpp
#define DOCINewtonOrbitalOptimizer_hpp

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "OrbitalOptimizationOptions.hpp"
#include "WaveFunction/WaveFunction.hpp"

#include <numopt.hpp>


namespace GQCP {


/**
 *  A class that performs gradient-and-Hessian-based orbital optimization for DOCI by sequentially
 *      - solving the DOCI eigenvalue problem
 *      - solving the Newton step to find the anti-Hermitian orbital rotation parameters
 *      - rotating the underlying spatial orbital basis
 */
class DOCINewtonOrbitalOptimizer {
private:
    GQCP::DOCI doci;  // the DOCI Hamiltonian builder
    GQCP::HamiltonianParameters ham_par;

    bool is_converged = false;
    std::vector<numopt::eigenproblem::Eigenpair> eigenpairs;  // eigenvalues and -vectors


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param doci instance and Hamiltonian parameters @param ham_par
     */
    DOCINewtonOrbitalOptimizer(const GQCP::DOCI& doci, const GQCP::HamiltonianParameters& ham_par);

    // GETTERS
    std::vector<numopt::eigenproblem::Eigenpair> get_eigenpairs() const;
    numopt::eigenproblem::Eigenpair get_eigenpair(size_t index = 0) const;

    // PUBLIC METHODS
    /**
     *  Perform the orbital optimization, given @param solver_options for the CI solver and the @param oo_options for the orbital optimization
     *
     *  The default values for the OrbitalOptimiationOptions are used when no options are supplied.
     */
    void solve(numopt::eigenproblem::BaseSolverOptions& solver_options, const GQCP::OrbitalOptimizationOptions& oo_options=GQCP::OrbitalOptimizationOptions());

    /**
     *  @return a WaveFunction instance after performing the orbital optimization for a given eigenvector at @param index
     */
    GQCP::WaveFunction get_wavefunction(size_t index = 0);
};


}  // namespace GQCP


#endif /* DOCINewtonOrbitalOptimizer_hpp */
