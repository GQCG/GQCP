#define BOOST_TEST_MODULE "CCD0"


#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/SpinorBasis/USpinorBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Mathematical/Representation/BlockRankFourTensor.hpp"
#include "Mathematical/Representation/QCRankFourTensor.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "QCMethod/CC/CCD0.hpp"
#include "QCMethod/HF/RHFSCFEnvironment.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/QCStructure.hpp"

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <fstream>
#include <string>



GQCP::ScalarSQTwoElectronOperator<double> twoElectronAntiSymmetrized(const GQCP::ScalarSQTwoElectronOperator<double>& g_op) {
    const auto& g_par = g_op.parameters();
    const size_t N = g_par.dimension(); // number of electrons
    
    GQCP::QCRankFourTensor<double> l_par(N);
    for (size_t p=0 ; p!=N ; p++){
        for (size_t q=0 ; q!=N ; q++){
            for (size_t r=0 ; r!=N ; r++){
                for (size_t s=0 ; s!=N ; s++){
                    l_par(p, q, r, s) = 2*g_par(p, q, r, s) - g_par(p, s, r, q);
                }
            }
        }
    }

    GQCP::ScalarSQTwoElectronOperator<double> l_op(l_par);

    return l_op;
}


double calculateProjection(const GQCP::ScalarSQOneElectronOperator<double>& f_op, const GQCP::ScalarSQTwoElectronOperator<double>& l_op, const GQCP::BlockRankFourTensor<double>& T, 
                            const size_t N, const size_t a, const size_t b, const size_t i, const size_t j){
    const auto& f_par = f_op.parameters();
    const auto& l_par = l_op.parameters();

    double projection = 0.0;
    // [ab||ij]
    projection += l_par(a/2, i/2, b/2, j/2);
    
    // t_ij^ac
    for (size_t c=N ; c!=2*N ; c++){
        projection += f_par(b/2, c/2) * T(a, c, i, j) - f_par(a/2, c/2) * T(b, c, i, j); // P(ab)
    }
    
    // t_ik^ab
    for (size_t k=0 ; k!=N ; k++){
        projection -= f_par(k/2, j/2) * T(a, b, i, k) - f_par(k/2, i/2) * T(a, b, j, k); // P(ij)
    }
    
    // t_ij^cd
    for (size_t c=N ; c!=2*N ; c++){ // electron 1
        for (size_t d=N ; d!=2*N ; d++){ // electron 2
            projection += 0.5 * l_par(a/2, c/2, b/2, d/2) * T(c, d, i, j);
        }
    }
    
    // t_kl^ab
    for (size_t k=0 ; k!=N ; k++){
        for (size_t l=0 ; l!=N ; l++){
            projection += 0.5 * l_par(k/2, i/2, l/2, j/2) * T(a, b, k, l);
        }
    }
    
    // t_ik^ac
    for (size_t k=0 ; k!=N ; k++){
        for (size_t c=N ; c!=2*N ; c++){
            projection += l_par(k/2, c/2, b/2, j/2) * T(a, c, i, k) // P(ab)
                        - l_par(k/2, c/2, a/2, j/2) * T(b, c, i, k)
                        - l_par(k/2, c/2, b/2, i/2) * T(a, c, j, k) // P(ij)
                        + l_par(k/2, c/2, a/2, i/2) * T(b, c, j, k);
        }
    }   
    
    // quadratic terms
    for (size_t k=0 ; k!=N ; k++){
        for (size_t c=N ; c!=2*N ; c++){
            for (size_t l=0 ; l!=N ; l++){
                for (size_t d=N ; d!=2*N ; d++){
                    // t_ik^ac * t_lj^db
                    projection += 0.5 * l_par(k/2, c/2, l/2, d/2) * T(a, c, i, k) * T(d, b, l, j) // P(ab)
                                - 0.5 * l_par(k/2, c/2, l/2, d/2) * T(b, c, i, k) * T(d, a, l, j) 
                                - 0.5 * l_par(k/2, c/2, l/2, d/2) * T(a, c, j, k) * T(d, b, l, i) // P(ij)
                                + 0.5 * l_par(k/2, c/2, l/2, d/2) * T(b, c, j, k) * T(d, a, l, i);

                    // t_ij^ac * t_kl^bd
                    projection -=  0.5 * l_par(k/2, c/2, l/2, d/2) * T(a, c, i, j) * T(b, d, k, l) // P(ab)
                                - 0.5 * l_par(k/2, c/2, l/2, d/2) * T(b, c, i, j) * T(a, d, k, l);

                    // t_ik^ab * t_jl^cd
                    projection -= 0.5 * l_par(k/2, c/2, l/2, d/2) * T(a, b, i, k) * T(c, d, j, l) // P(ij)
                                - 0.5 * l_par(k/2, c/2, l/2, d/2) * T(a, b, j, k) * T(c, d, i, l);

                    // t_ij^cd * t_kl^ab
                    projection += 0.25 * l_par(k/2, c/2, l/2, d/2) * T(c, d, i, j) * T(a, b, k, l);
                }
            }
        }
    }
    
    
    return projection;
}

GQCP::BlockRankFourTensor<double> calculateProjections(const GQCP::ScalarSQOneElectronOperator<double>& f_op, const GQCP::ScalarSQTwoElectronOperator<double>& l_op, 
                                                        const GQCP::BlockRankFourTensor<double>& T, const size_t N){
    size_t occupied_s = 0, occupied_e = N;
    size_t virtual_s = N, virtual_e = 2*N;
    GQCP::BlockRankFourTensor<double> F(virtual_s, virtual_e, virtual_s, virtual_e, occupied_s, occupied_e, occupied_s, occupied_e);

    for (size_t a=N ; a!=2*N ; a++){
        for (size_t b=N ; b!=2*N ; b++){
            for (size_t i=0 ; i!=N ; i++){
                for (size_t j=0 ; j!=N ; j++){
                    F(a, b, i, j) = calculateProjection(f_op, l_op, T, N, a, b, i, j);
                }
            }
        }
    }

    return F;
}

GQCP::BlockRankFourTensor<double> initializeAmplitudes(const GQCP::ScalarSQOneElectronOperator<double>& f_op, const GQCP::ScalarSQTwoElectronOperator<double>& l_op, const size_t N){

    const auto& f_par = f_op.parameters();
    const auto& l_par = l_op.parameters();

    size_t occupied_s = 0, occupied_e = N;
    size_t virtual_s = N, virtual_e = 2*N;
    GQCP::BlockRankFourTensor<double> T(virtual_s, virtual_e, virtual_s, virtual_e, occupied_s, occupied_e, occupied_s, occupied_e);

    // initialize t-values
    for (size_t a=N ; a!=2*N ; a++){
        for (size_t b=N ; b!=2*N ; b++){
            for (size_t i=0 ; i!=N ; i++){
                for (size_t j=0 ; j!=N ; j++){
                    //std::cout<<"<"<<a<<b<<"|"<<i<<j<<">"<<std::endl;
                    //std::cout<<"("<<a<<i<<"|"<<b<<j<<")"<<std::endl;
                    //std::cout<<ls(a/2, i/2, b/2, j/2)<<"/("<<fs(i/2, i/2)<<"+"<<fs(j/2, j/2)<<"-"<<fs(a/2, a/2)<<"-"<<fs(b/2, b/2)<<")"<<std::endl;
                    //std::cout<<ls(a/2, i/2, b/2, j/2) / (fs(i/2, i/2) + fs(j/2, j/2) - fs(a/2, a/2) - fs(b/2, b/2))<<std::endl;
                    T(a, b, i, j) = l_par(a/2, i/2, b/2, j/2) / (f_par(i/2, i/2) + f_par(j/2, j/2) - f_par(a/2, a/2) - f_par(b/2, b/2));
                }
            }
        }
    }
    return T;
}

void step(const GQCP::ScalarSQOneElectronOperator<double>& f_op, const GQCP::ScalarSQTwoElectronOperator<double>& l_op, 
            GQCP::BlockRankFourTensor<double>& T, const GQCP::BlockRankFourTensor<double>& F, const size_t N){
    
    const auto& f_par = f_op.parameters();

    for (size_t a=N ; a!=2*N ; a++){
        for (size_t b=N ; b!=2*N ; b++){
            for (size_t i=0 ; i!=N ; i++){
                for (size_t j=0 ; j!=N ; j++){
                    T(a, b, i, j) += F(a, b, i, j) / (f_par(i/2, i/2) + f_par(j/2, j/2) + f_par(a/2, a/2) + f_par(b/2, b/2));
                }
            }
        }
    }
}

double calculateEnergy(const double e_hf, const GQCP::ScalarSQTwoElectronOperator<double>& l_op, 
            GQCP::BlockRankFourTensor<double>& T, const GQCP::BlockRankFourTensor<double>& F, const size_t N){
    const auto& l_par = l_op.parameters();
    
    double energy = 0.0;

    for (size_t i=0 ; i!=N ; i++){
        for (size_t j=0 ; j!=N ; j++){
            for (size_t a=N ; a!=2*N ; a++){
                for (size_t b=N ; b!=2*N ; b++){
                    energy += 0.25 * l_par(i/2, a/2, j/2, b/2) * T(a, b, i, j);
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
    GQCP::ScalarSQOneElectronOperator<double> inactive_f_op = sq_hamiltonian.calculateInactiveFockian(N/2);
    
    // setting up the environment and solver for the RHF energy
    GQCP::RHFSCFEnvironment<double> environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, ss);
    GQCP::IterativeAlgorithm<GQCP::RHFSCFEnvironment<double>> solver = GQCP::RHFSCFSolver<double>::DIIS();
    GQCP::DiagonalRHFFockMatrixObjective<double> objective(sq_hamiltonian);
    
    // calculate necessary E_HF
    GQCP::QCMethod::RHF<double> rhf;
    GQCP::QCStructure<GQCP::QCModel::RHF<double>> rhf_parameters = rhf.optimize(objective, solver, environment);
    double e_rhf = rhf_parameters.groundStateEnergy();
    std::cout<<"RHF ground state energy: "<<e_rhf<<" hartree"<<std::endl;

    // calculate antisymmetrized two electron integrals
    const auto l_op = twoElectronAntiSymmetrized(sq_hamiltonian.twoElectron());

    // initialize t amplitudes
    auto T = initializeAmplitudes(inactive_f_op, l_op, N);
    T.asTensor().print();
    std::cout<<"\n\n";
    
    // calculate the function values for the initial t amplitudes
    auto F = calculateProjections(inactive_f_op, l_op, T, N);
    F.asTensor().print();
    
    // calculate the CCD energy
    auto E_old = calculateEnergy(e_rhf, l_op, T, F, N);
    std::cout<<E_old;
    
    for (int s=0 ; s<50 ; s++){
        step(inactive_f_op, l_op, T, F, N);
        F = calculateProjections(inactive_f_op, l_op, T, N);

        std::cout << calculateEnergy(e_rhf, l_op, T, F, N) << std::endl;
    }
    
 
}
