#define BOOST_TEST_MODULE "CCD0"

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>
#include <string>

#include "QCMethod/CC/CCD0.hpp"
#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/SpinorBasis/USpinorBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Mathematical/Representation/QCRankFourTensor.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "QCMethod/HF/RHFSCFEnvironment.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/QCStructure.hpp"
#include "Mathematical/Representation/BlockRankFourTensor.hpp"

GQCP::QCRankFourTensor<double> twoElectronAntiSymmetrized(const GQCP::QCRankFourTensor<double>& gs) {
    size_t dim = gs.dimension();
    GQCP::QCRankFourTensor<double> ls(dim);

    for (size_t p=0 ; p!=dim ; p++){
        for (size_t q=0 ; q!=dim ; q++){
            for (size_t r=0 ; r!=dim ; r++){
                for (size_t s=0 ; s!=dim ; s++){
                    ls(p, q, r, s) = 2*gs(p, q, r, s) - gs(p, s, r, q);
                }
            }
        }
    }
    return ls;
}

double calculateProjection(const GQCP::SQHamiltonian<double>& sq_hamiltonian, GQCP::BlockRankFourTensor<double>& ts, const GQCP::ScalarSQOneElectronOperator<double>& f,
                        const size_t a, const size_t b, const size_t i, const size_t j){
    const auto hs = sq_hamiltonian.core().parameters(0);
    const auto gs = sq_hamiltonian.twoElectron().parameters(0);
    const auto ls = twoElectronAntiSymmetrized(gs);
    const auto fs = f.parameters();

    const size_t N = gs.dimension(); // returns the number of electrons

    double projection = 0.0;
    // [ab||ij]
    projection += ls(a/2, i/2, b/2, j/2);
    /*
    // t_ij^ac
    for (size_t c=N ; c!=2*N ; c++){
        projection += fs(b/2, c/2) * ts(a, c, i, j) - fs(a/2, c/2) * ts(b, c, i, j); // P(ab)
    }

    // t_ik^ab
    for (size_t k=0 ; k!=N ; k++){
        projection -= fs(k/2, j/2) * ts(a, b, i, k) - fs(k/2, i/2) * ts(a, b, j, k); // P(ij)
    }
    
    // t_ij^cd
    for (size_t c=N ; c!=2*N ; c++){ // electron 1
        for (size_t d=N ; d!=2*N ; d++){ // electron 2
            projection += 0.5 * ls(a/2, b/2, c/2, d/2) * ts(c, d, i, j);
        }
    }

    // t_kl^ab
    for (size_t k=0 ; k!=N ; k++){
        for (size_t l=0 ; l!=N ; l++){
            projection += 0.5 * ls(k/2, l/2, i/2, j/2) * ts(a, b, k, l);
        }
    }

    // t_ik^ac
    for (size_t k=0 ; k!=N ; k++){
        for (size_t c=N ; c!=2*N ; c++){
            projection += ls(k/2, b/2, c/2, j/2) * ts(a, c, i, k) // P(ab)
                        - ls(k/2, a/2, c/2, j/2) * ts(b, c, i, k)
                        - ls(k/2, b/2, c/2, i/2) * ts(a, c, j, k) // P(ij)
                        + ls(k/2, a/2, c/2, i/2) * ts(b, c, j, k);
        }
    }   
    
    // quadratic terms
    for (size_t k=0 ; k!=N ; k++){
        for (size_t c=N ; c!=2*N ; c++){
            for (size_t l=0 ; l!=N ; l++){
                for (size_t d=N ; d!=2*N ; d++){
                    // t_ik^ac * t_lj^db
                    projection += 0.5 * ls(k/2, l/2, c/2, d/2) * ts(a, c, i, k) * ts(d, b, l, j) // P(ab)
                                - 0.5 * ls(k/2, l/2, c/2, d/2) * ts(b, c, i, k) * ts(d, a, l, j) 
                                - 0.5 * ls(k/2, l/2, c/2, d/2) * ts(a, c, j, k) * ts(d, b, l, i) // P(ij)
                                + 0.5 * ls(k/2, l/2, c/2, d/2) * ts(b, c, j, k) * ts(d, a, l, i);

                    // t_ik^cd * t_lj^ab
                    projection +=  0.5 * ls(k/2, l/2, c/2, d/2) * ts(c, d, i, k) * ts(a, b, l, j) // P(ij)
                                - 0.5 * ls(k/2, l/2, c/2, d/2) * ts(c, d, j, k) * ts(a, b, l, i);

                    // t_kl^ac * t_ij^db
                    projection += 0.5 * ls(k/2, l/2, c/2, d/2) * ts(a, c, k, l) * ts(d, b, i, j) // P(ab)
                                - 0.5 * ls(k/2, l/2, c/2, d/2) * ts(b, c, k, l) * ts(d, a, i, j);

                    // t_ij^cd * t_kl^ab
                    projection += 0.25 * ls(k/2, l/2, c/2, d/2) * ts(c, d, i, j) * ts(a, b, k, l);
                }
            }
        }
    }
    */
    
    return projection;
}

GQCP::BlockRankFourTensor<double> calculateProjections(const GQCP::SQHamiltonian<double>& sq_hamiltonian, GQCP::BlockRankFourTensor<double>& ts, const GQCP::ScalarSQOneElectronOperator<double>& f){
    const size_t N = sq_hamiltonian.core().parameters(0).dimension(); // returns the number of electrons

    size_t occupied_s = 0, occupied_e = N;
    size_t virtual_s = N, virtual_e = 2*N;

    GQCP::BlockRankFourTensor<double> F(virtual_s, virtual_e, virtual_s, virtual_e, occupied_s, occupied_e, occupied_s, occupied_e);

    for (size_t a=N ; a!=2*N ; a++){
        for (size_t b=N ; b!=2*N ; b++){
            for (size_t i=0 ; i!=N ; i++){
                for (size_t j=0 ; j!=N ; j++){
                    F(a, b, i, j) = calculateProjection(sq_hamiltonian, ts, f, a, b, i, j);
                }
            }
        }
    }

    return F;
}

GQCP::BlockRankFourTensor<double> initializeAmplitudes(const GQCP::SQHamiltonian<double>& sq_hamiltonian, const GQCP::ScalarSQOneElectronOperator<double>& f){
    const auto gs = sq_hamiltonian.twoElectron().parameters(0);
    const auto ls = twoElectronAntiSymmetrized(gs);

    const size_t N = gs.dimension(); // number of electrons
    const size_t K = N*2; // number of spinors
    const auto fs = f.parameters();
    //const auto hs = sq_hamiltonian.core().parameters(0);

    size_t occupied_s = 0, occupied_e = N;
    size_t virtual_s = N, virtual_e = 2*N;
    GQCP::BlockRankFourTensor<double> ts(virtual_s, virtual_e, virtual_s, virtual_e, occupied_s, occupied_e, occupied_s, occupied_e);

    // initialize t-values
    for (size_t a=N ; a!=2*N ; a++){
        for (size_t b=N ; b!=2*N ; b++){
            for (size_t i=0 ; i!=N ; i++){
                for (size_t j=0 ; j!=N ; j++){
                    //std::cout<<"<"<<a<<b<<"|"<<i<<j<<">"<<std::endl;
                    //std::cout<<"("<<a<<i<<"|"<<b<<j<<")"<<std::endl;
                    //std::cout<<ls(a/2, i/2, b/2, j/2)<<"/("<<fs(i/2, i/2)<<"+"<<fs(j/2, j/2)<<"-"<<fs(a/2, a/2)<<"-"<<fs(b/2, b/2)<<")"<<std::endl;
                    //std::cout<<ls(a/2, i/2, b/2, j/2) / (fs(i/2, i/2) + fs(j/2, j/2) - fs(a/2, a/2) - fs(b/2, b/2))<<std::endl;
                    ts(a, b, i, j) = ls(a/2, i/2, b/2, j/2) / (fs(i/2, i/2) + fs(j/2, j/2) - fs(a/2, a/2) - fs(b/2, b/2));
                }
            }
        }
    }
    return ts;
}

void step(const GQCP::SQHamiltonian<double>& sq_hamiltonian, GQCP::BlockRankFourTensor<double>& ts, GQCP::BlockRankFourTensor<double>& F, const GQCP::ScalarSQOneElectronOperator<double>& f){
    const auto fs = f.parameters();
    const size_t N = sq_hamiltonian.core().parameters(0).dimension(); // returns the number of electrons

    for (size_t a=N ; a!=2*N ; a++){
        for (size_t b=N ; b!=2*N ; b++){
            for (size_t i=0 ; i!=N ; i++){
                for (size_t j=0 ; j!=N ; j++){
                    ts(a, b, i, j) += F(a, b, i, j) / (fs(i/2, i/2) + fs(j/2, j/2) + fs(a/2, a/2) + fs(b/2, b/2));
                }
            }
        }
    }
}

double calculateEnergy(const double e_hf, const GQCP::SQHamiltonian<double>& sq_hamiltonian, GQCP::BlockRankFourTensor<double>& ts){
    const auto hs = sq_hamiltonian.core().parameters(0);
    const auto gs = sq_hamiltonian.twoElectron().parameters(0);
    const auto ls = gs; //twoElectronAntiSymmetrized(gs);

    const size_t N = hs.dimension(); // returns the number of electrons
    
    double energy = 0.0;

    for (size_t i=0 ; i!=N ; i++){
        for (size_t j=0 ; j!=N ; j++){
            for (size_t a=N ; a!=2*N ; a++){
                for (size_t b=N ; b!=2*N ; b++){
                    energy += 0.25 * ls(i/2, a/2, j/2, b/2) * ts(a, b, i, j);
                }
            }
        }
    }
    
    energy += e_hf;
    
    return energy;
}

BOOST_AUTO_TEST_CASE ( CCD0 ) {

    const auto molecule = GQCP::Molecule::ReadXYZ("/Users/daria/ugent/GQCP/gqcp/tests/data/h2.xyz" , 0);
    const size_t N = molecule.numberOfElectrons();

    GQCP::RSpinorBasis<double, GQCP::GTOShell> rspinor_basis(molecule, "STO-3G");
    std::cout<<"Number of spatial orbitals: "<<rspinor_basis.numberOfSpatialOrbitals()<<std::endl;
    std::cout<<"Number of spinors: "<<rspinor_basis.numberOfSpinors()<<std::endl;


    // calculate necessary integrals and operators
    const auto S = rspinor_basis.quantize(GQCP::Operator::Overlap());
    const auto ss = S.parameters(0);
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(rspinor_basis, molecule);
    GQCP::ScalarSQOneElectronOperator<double> inactive_fock_matrix = sq_hamiltonian.calculateInactiveFockian(N/2);
    
    // setting up the environment and solver for the RHF energy
    GQCP::RHFSCFEnvironment<double> environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, ss);
    GQCP::IterativeAlgorithm<GQCP::RHFSCFEnvironment<double>> solver = GQCP::RHFSCFSolver<double>::DIIS();
    GQCP::DiagonalRHFFockMatrixObjective<double> objective(sq_hamiltonian);
    
    // calculate necessary E_HF
    GQCP::QCMethod::RHF<double> rhf;
    GQCP::QCStructure<GQCP::QCModel::RHF<double>> rhf_parameters = rhf.optimize(objective, solver, environment);
    double rhf_energy = rhf_parameters.groundStateEnergy();
    std::cout<<"RHF ground state energy: "<<rhf_energy<<" hartree"<<std::endl;
    
    // initialize t amplitudes
    auto ts = initializeAmplitudes(sq_hamiltonian, inactive_fock_matrix);
    ts.asTensor().print();
    std::cout<<"\n\n";

    // calculate the function values for the initial t amplitudes
    auto F = calculateProjections(sq_hamiltonian, ts, inactive_fock_matrix);
    F.asTensor().print();

    //for (int s=0 ; s<10 ; s++){
    //    step(sq_hamiltonian, ts, F, inactive_fock_matrix);
    //    F = calculateProjections(sq_hamiltonian, ts, inactive_fock_matrix);
    //    std::cout<<calculateEnergy(rhf_energy, sq_hamiltonian, ts)<<std::endl;
    //}
 
}
