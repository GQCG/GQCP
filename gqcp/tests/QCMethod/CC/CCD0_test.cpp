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
#include <stdlib.h>



GQCP::ScalarSQTwoElectronOperator<double> twoElectronAntiSymmetrized(const GQCP::ScalarSQTwoElectronOperator<double>& g_op) {
    const auto& g_par = g_op.parameters();
    const size_t K = g_par.dimension(); // number of spatial orbitals
    
    GQCP::QCRankFourTensor<double> l_par(K);
    for (size_t p=0 ; p!=K ; p++){
        for (size_t q=0 ; q!=K ; q++){
            for (size_t r=0 ; r!=K ; r++){
                for (size_t s=0 ; s!=K ; s++){
                    l_par(p, q, r, s) = 2*g_par(p, q, r, s) - g_par(p, s, r, q);
                }
            }
        }
    }

    GQCP::ScalarSQTwoElectronOperator<double> l_op(l_par);

    return l_op;
}


double calculateProjection(const GQCP::ScalarSQOneElectronOperator<double>& f_op, const GQCP::ScalarSQTwoElectronOperator<double>& l_op, const GQCP::BlockRankFourTensor<double>& T, 
                            const size_t N_P, const size_t a, const size_t b, const size_t i, const size_t j){
    const auto& f_par = f_op.parameters();
    const auto& l_par = l_op.parameters();

    const size_t K = f_par.dimension();

    double projection = 0.0;
    // [ab||ij]
    projection += l_par(a, i, b, j);
    
    // t_ij^ac
    for (size_t c=N_P ; c!=K ; c++){
        projection += f_par(b, c) * T(a, c, i, j) - f_par(a, c) * T(b, c, i, j); // P(ab)
    }
    
    // t_ik^ab
    for (size_t k=0 ; k!=N_P ; k++){
        projection -= f_par(k, j) * T(a, b, i, k) - f_par(k, i) * T(a, b, j, k); // P(ij)
    }

    // t_kl^ab
    for (size_t k=0 ; k!=N_P ; k++){
        for (size_t l=0 ; l!=N_P ; l++){
            projection += 0.5 * l_par(k, i, l, j) * T(a, b, k, l);
        }
    }
    
    // t_ij^cd
    for (size_t c=N_P ; c!=K ; c++){ // electron 1
        for (size_t d=N_P ; d!=K ; d++){ // electron 2
            projection += 0.5 * l_par(a, c, b, d) * T(c, d, i, j);
        }
    }
    
    // t_ik^ac
    for (size_t k=0 ; k!=N_P ; k++){
        for (size_t c=N_P ; c!=K ; c++){
            projection += l_par(k, c, b, j) * T(a, c, i, k) // P(ab)
                        - l_par(k, c, a, j) * T(b, c, i, k)
                        - l_par(k, c, b, i) * T(a, c, j, k) // P(ij)
                        + l_par(k, c, a, i) * T(b, c, j, k);
        }
    }   
    
    // quadratic terms
    for (size_t k=0 ; k!=N_P ; k++){
        for (size_t c=N_P ; c!=K ; c++){
            for (size_t l=0 ; l!=N_P ; l++){
                for (size_t d=N_P ; d!=K ; d++){
                    // t_ik^ac * t_lj^db
                    projection += 0.5 * l_par(k, c, l, d) * T(a, c, i, k) * T(d, b, l, j) // P(ab)
                                - 0.5 * l_par(k, c, l, d) * T(b, c, i, k) * T(d, a, l, j) 
                                - 0.5 * l_par(k, c, l, d) * T(a, c, j, k) * T(d, b, l, i) // P(ij)
                                + 0.5 * l_par(k, c, l, d) * T(b, c, j, k) * T(d, a, l, i);

                    // t_ij^cd * t_kl^ab
                    projection += 0.25 * l_par(k, c, l, d) * T(c, d, i, j) * T(a, b, k, l);

                    // t_ij^ac * t_kl^bd
                    projection -=  0.5 * l_par(k, c, l, d) * T(a, c, i, j) * T(b, d, k, l) // P(ab)
                                - 0.5 * l_par(k, c, l, d) * T(b, c, i, j) * T(a, d, k, l);

                    // t_ik^ab * t_jl^cd
                    projection -= 0.5 * l_par(k, c, l, d) * T(a, b, i, k) * T(c, d, j, l) // P(ij)
                                - 0.5 * l_par(k, c, l, d) * T(a, b, j, k) * T(c, d, i, l);
                }
            }
        }
    }
    
    return projection;
}

GQCP::BlockRankFourTensor<double> calculateProjections(const GQCP::ScalarSQOneElectronOperator<double>& f_op, const GQCP::ScalarSQTwoElectronOperator<double>& l_op, 
                                                        const GQCP::BlockRankFourTensor<double>& T, const size_t N_P){
    const size_t K = f_op.parameters().dimension();

    GQCP::BlockRankFourTensor<double> F(N_P, K, N_P, K, 0, N_P, 0, N_P); // virtual-virtual-occupied-occupied

    for (size_t a=N_P ; a!=K ; a++){
        for (size_t b=N_P ; b!=K ; b++){
            for (size_t i=0 ; i!=N_P ; i++){
                for (size_t j=0 ; j!=N_P ; j++){
                    F(a, b, i, j) = calculateProjection(f_op, l_op, T, N_P, a, b, i, j);
                }
            }
        }
    }

    return F;
}

GQCP::BlockRankFourTensor<double> initializeAmplitudes(const GQCP::ScalarSQOneElectronOperator<double>& f_op, const GQCP::ScalarSQTwoElectronOperator<double>& l_op, const size_t N_P){

    const auto& f_par = f_op.parameters();
    const auto& l_par = l_op.parameters();

    const size_t K = f_par.dimension();
    //std::cout<<"N_P: "<<N_P<<"  K: "<<K<<std::endl;

    GQCP::BlockRankFourTensor<double> T(N_P, K, N_P, K, 0, N_P, 0, N_P);

    // initialize t-values
    for (size_t a=N_P ; a!=K ; a++){
        for (size_t b=N_P ; b!=K ; b++){
            for (size_t i=0 ; i!=N_P ; i++){
                for (size_t j=0 ; j!=N_P ; j++){
                    //std::cout<<"<"<<a<<b<<"|"<<i<<j<<">"<<std::endl;
                    //std::cout<<"("<<a<<i<<"|"<<b<<j<<")"<<std::endl;
                    //std::cout<<ls(a/2, i/2, b/2, j/2)<<"/("<<fs(i/2, i/2)<<"+"<<fs(j/2, j/2)<<"-"<<fs(a/2, a/2)<<"-"<<fs(b/2, b/2)<<")"<<std::endl;
                    //std::cout<<ls(a/2, i/2, b/2, j/2) / (fs(i/2, i/2) + fs(j/2, j/2) - fs(a/2, a/2) - fs(b/2, b/2))<<std::endl;
                    T(a, b, i, j) = l_par(a, i, b, j) / (f_par(i, i) + f_par(j, j) - f_par(a, a) - f_par(b, b));
                }
            }
        }
    }
    return T;
}

void step(const GQCP::ScalarSQOneElectronOperator<double>& f_op, const GQCP::ScalarSQTwoElectronOperator<double>& l_op, 
            GQCP::BlockRankFourTensor<double>& T, const GQCP::BlockRankFourTensor<double>& F, const size_t N_P){
    
    const auto& f_par = f_op.parameters();

    const size_t K = f_par.dimension();

    for (size_t a=N_P ; a!=K ; a++){
        for (size_t b=N_P ; b!=K ; b++){
            for (size_t i=0 ; i!=N_P ; i++){
                for (size_t j=0 ; j!=N_P ; j++){
                    //std::cout<<"T("<<a<<b<<i<<j<<") = "<<F(a, b, i, j)<<" / "<<f_par(i, i) + f_par(j, j) + f_par(a, a) + f_par(b, b)<<std::endl;
                    T(a, b, i, j) += F(a, b, i, j) / (f_par(i, i) + f_par(j, j) + f_par(a, a) + f_par(b, b));
                }
            }
        }
    }
}

double calculateEnergy(const double e_hf, const GQCP::ScalarSQTwoElectronOperator<double>& l_op, 
            GQCP::BlockRankFourTensor<double>& T, const GQCP::BlockRankFourTensor<double>& F, const size_t N_P){
    const auto& l_par = l_op.parameters();
    
    const size_t K = l_par.dimension();

    double energy = 0.0;

    for (size_t a=N_P ; a!=K ; a++){
        for (size_t b=N_P ; b!=K ; b++){
            for (size_t i=0 ; i!=N_P ; i++){
                for (size_t j=0 ; j!=N_P ; j++){
                    energy += 0.25 * l_par(i, a, j, b) * T(a, b, i, j);
                }
            }
        }
    }
    
    //energy += e_hf;
    
    return energy;
}

BOOST_AUTO_TEST_CASE ( CCD0 ) {

    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2.xyz" , 0);
    const size_t N = molecule.numberOfElectrons();
    const size_t N_P = N/2;

    GQCP::RSpinorBasis<double, GQCP::GTOShell> rspinor_basis(molecule, "STO-3G");
    std::cout<<"Number of spatial orbitals: "<<rspinor_basis.numberOfSpatialOrbitals()<<std::endl;
    std::cout<<"Number of spinors: "<<rspinor_basis.numberOfSpinors()<<std::endl;


    // calculate necessary integrals and operators
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(rspinor_basis, molecule);
    GQCP::ScalarSQOneElectronOperator<double> inactive_f_op = sq_hamiltonian.calculateInactiveFockian(N_P);
    
    // setting up the environment and solver for the RHF energy
    GQCP::RHFSCFEnvironment<double> environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, rspinor_basis.overlap().parameters());
    GQCP::IterativeAlgorithm<GQCP::RHFSCFEnvironment<double>> solver = GQCP::RHFSCFSolver<double>::Plain();
    //GQCP::DiagonalRHFFockMatrixObjective<double> objective(sq_hamiltonian);
    solver.perform(environment);

    // check the total energy
    const double e_rhf = environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(molecule).value();
    // calculate necessary E_HF
    //GQCP::QCMethod::RHF<double> rhf;
    //GQCP::QCStructure<GQCP::QCModel::RHF<double>> rhf_parameters = rhf.optimize(objective, solver, environment);
    //double e_rhf = rhf_parameters.groundStateEnergy();
    std::cout<<"RHF ground state energy: "<<e_rhf<<" hartree.\n\n"<<std::endl;

    // calculate antisymmetrized two electron integrals
    const auto l_op = twoElectronAntiSymmetrized(sq_hamiltonian.twoElectron());
    
    // initialize t amplitudes
    auto T = initializeAmplitudes(inactive_f_op, l_op, N_P);
    
    // calculate the function values for the initial t amplitudes
    auto F = calculateProjections(inactive_f_op, l_op, T, N_P);
    
    // calculate the CCD energy
    double E_new;
    auto E_old = calculateEnergy(e_rhf, l_op, T, F, N_P);
    std::cout<<"E (initial): "<<E_old<<std::endl;
    
    for (int s=0 ; s<100 ; s++){
        step(inactive_f_op, l_op, T, F, N_P);
        F = calculateProjections(inactive_f_op, l_op, T, N_P);
        E_new = calculateEnergy(e_rhf, l_op, T, F, N_P);
        //std::cout << E_new << std::endl;

        if (abs(E_new-E_old)<10e-4){
            std::cout<<"Convergence after "<<s<<" iterations."<<std::endl;
            break;
        }
    }
    std::cout<<"The final, converged energy is: "<<E_new<<std::endl;
    
 
}
