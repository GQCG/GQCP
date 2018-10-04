#ifndef PlainRHFSCFSolver_hpp
#define PlainRHFSCFSolver_hpp


#include "RHFSCFSolver.hpp"


namespace GQCG {


/**
 *  A plain RHF SCF solver
 */
class PlainRHFSCFSolver : public GQCG::RHFSCFSolver {
private:
    /**
     *  Calculate a new Fock matrix (in AO basis)
     */
    Eigen::MatrixXd calculateNewFockMatrix(const Eigen::MatrixXd& D_AO) override;

public:
    // CONSTRUCTORS
    /**
     *  Constructor based on given @param hamiltonian parameters, @param molecule, @param maximum_number_of_iterations and @param SCF threshold
     */
    PlainRHFSCFSolver(GQCG::HamiltonianParameters ham_par, GQCG::Molecule molecule, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);

};

}  // namespace GQCG


#endif /* PlainRHFSCFSolver_hpp */
