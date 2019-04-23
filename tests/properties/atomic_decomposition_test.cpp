// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2019  the GQCG developers
//
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
//
#define BOOST_TEST_MODULE "atomic_decomposition"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "RHF/PlainRHFSCFSolver.hpp"
#include "CISolver/CISolver.hpp"
#include "RDM/RDMCalculator.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "properties/expectation_values.hpp"
#include "units.hpp"

GQCP::TwoElectronOperator<double> transform_rdm(const GQCP::TwoRDM<double>& x, const GQCP::SquareMatrix<double>& T) {

    // Since we're only getting T as a matrix, we should make the appropriate tensor to perform contractions
    // For the const argument, we need the const in the template
    //      For more info, see: https://stackoverflow.com/questions/45283468/eigen-const-tensormap
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> T_tensor (T.data(), T.rows(), T.cols());


    // We will have to do four single contractions, so we specify the contraction indices
    // Eigen3 does not document its tensor contraction clearly, so see the accepted answer on stackoverflow (https://stackoverflow.com/a/47558349/7930415):
    //      Eigen3 does not accept a way to specify the output axes: instead, it retains the order from left to right of the axes that survive the contraction.
    //      This means that, in order to get the right ordering of the axes, we will have to swap axes

    // g(T U V W)  T^*(V R) -> a(T U R W) but we get a(T U W R)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair1 = {Eigen::IndexPair<int>(2, 0)};
    Eigen::array<int, 4> shuffle_1 {0, 1, 3, 2};

    // a(T U R W)  T(W S) -> b(T U R S) and we get b(T U R S), so no shuffle is needed
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair2 = {Eigen::IndexPair<int>(3, 0)};

    // T(U Q)  b(T U R S) -> c(T Q R S) but we get c(Q T R S)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair3 = {Eigen::IndexPair<int>(0, 1)};
    Eigen::array<int, 4> shuffle_3 {1, 0, 2, 3};

    // T^*(T P)  c(T Q R S) -> g'(P Q R S) and we get g_SO(P Q R S), so no shuffle is needed
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair4 = {Eigen::IndexPair<int>(0, 0)};


    // Calculate the contractions. We write this as one large contraction to
    //  1) avoid storing intermediate contractions
    //  2) let Eigen figure out some optimizations
    GQCP::TwoElectronOperator<double> g_transformed = T_tensor.conjugate().contract(T_tensor.contract(x.contract(T_tensor.conjugate(), contraction_pair1).shuffle(shuffle_1).contract(T_tensor, contraction_pair2), contraction_pair3).shuffle(shuffle_3), contraction_pair4);

    return g_transformed;
}





GQCP::TwoElectronOperator<double> transform_rdm(const GQCP::TwoRDM<double>& x, const GQCP::SquareMatrix<double>& ls, const GQCP::SquareMatrix<double>& rs) {

    // Since we're only getting T as a matrix, we should make the appropriate tensor to perform contractions
    // For the const argument, we need the const in the template
    //      For more info, see: https://stackoverflow.com/questions/45283468/eigen-const-tensormap

    Eigen::MatrixXd id = Eigen::MatrixXd::Identity(ls.rows(), ls.cols());
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> ls_tensor (ls.data(), ls.rows(), ls.cols());
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> rs_tensor (rs.data(), ls.rows(), ls.cols());
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> identity (id.data(), ls.rows(), ls.cols());


    // We will have to do four single contractions, so we specify the contraction indices
    // Eigen3 does not document its tensor contraction clearly, so see the accepted answer on stackoverflow (https://stackoverflow.com/a/47558349/7930415):
    //      Eigen3 does not accept a way to specify the output axes: instead, it retains the order from left to right of the axes that survive the contraction.
    //      This means that, in order to get the right ordering of the axes, we will have to swap axes

    // g(T U V W)  T^*(V R) -> a(T U R W) but we get a(T U W R)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair1 = {Eigen::IndexPair<int>(2, 0)};
    Eigen::array<int, 4> shuffle_1 {0, 1, 3, 2};

    // a(T U R W)  T(W S) -> b(T U R S) and we get b(T U R S), so no shuffle is needed
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair2 = {Eigen::IndexPair<int>(3, 0)};

    // T(U Q)  b(T U R S) -> c(T Q R S) but we get c(Q T R S)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair3 = {Eigen::IndexPair<int>(0, 1)};
    Eigen::array<int, 4> shuffle_3 {1, 0, 2, 3};

    // T^*(T P)  c(T Q R S) -> g'(P Q R S) and we get g_SO(P Q R S), so no shuffle is needed
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair4 = {Eigen::IndexPair<int>(0, 0)};


    // Calculate the contractions. We write this as one large contraction to
    //  1) avoid storing intermediate contractions
    //  2) let Eigen figure out some optimizations
    GQCP::TwoElectronOperator<double> g_transformed = ls_tensor.conjugate().
            contract(identity.
            contract(x.contract(rs_tensor.conjugate(), contraction_pair1).shuffle(shuffle_1).
            contract(identity, contraction_pair2), contraction_pair3).shuffle(shuffle_3), contraction_pair4);;

    return g_transformed;
}


BOOST_AUTO_TEST_CASE ( decomposition_CO_STO_3G ) {

    // Create the molecular Hamiltonian parameters in an AO basis
    //auto CO = GQCP::Molecule::Readxyz("data/CO_mulliken.xyz");
    auto CO = GQCP::Molecule::Readxyz("data/h2o.xyz");
    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(CO, "STO-3G");
    auto ham_par_ao_copy = mol_ham_par;
    auto K = mol_ham_par.get_K();

    GQCP::Vectoru gto_list (K-2);
    GQCP::Vectoru gto_list3 (2);

    Eigen::MatrixXd pa = Eigen::MatrixXd::Zero(K,K);
    Eigen::MatrixXd pb = Eigen::MatrixXd::Zero(K,K);
    for(size_t i = 0; i<K-2; i++){
        gto_list[i] = i;
        std::cout<<i<<" ";
        pa(i,i) = 1;
    }
    std::cout<<std::setprecision(16);
    for(size_t i = K-2; i<K; i++){
        gto_list3[i-K+2] = i;
        std::cout<<i<<" ";
        pb(i,i) = 1;
    }
    std::cout<<std::endl;

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, CO);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    const auto& T = rhf.get_C();

    // Transform the ham_par
    mol_ham_par.transform(T);

    GQCP::ProductFockSpace fock_space (K, CO.get_N()/2, CO.get_N()/2);  // dim = 441

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, mol_ham_par);

    // Solve Davidson
    GQCP::DavidsonSolverOptions solver_options(fock_space.HartreeFockExpansion());
    ci_solver.solve(solver_options);

    // Retrieve the eigenvector
    auto fci_coeff = ci_solver.get_eigenpair().get_eigenvector();
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    GQCP::RDMCalculator rdm_calc(fock_space);
    rdm_calc.set_coefficients(fci_coeff);

    auto one_rdm = rdm_calc.calculate1RDMs().one_rdm;
    auto two_rdm = rdm_calc.calculate2RDMs().two_rdm;

    // Convert rdm to AObasis
    GQCP::OneRDM<double> ao_one_rdm = T * (one_rdm) * T.adjoint() ;
    GQCP::TwoRDM<double> ao_two_rdm = transform_rdm(two_rdm, T.adjoint());

    double trace_energy = GQCP::calculateExpectationValue(ham_par_ao_copy, ao_one_rdm, ao_two_rdm);
    double lol = 0;


    auto S = ham_par_ao_copy.get_S();
    auto h = ham_par_ao_copy.get_h();
    auto h2 = mol_ham_par.get_h();

    auto g = ham_par_ao_copy.get_g();
    auto g2 = mol_ham_par.get_g();

    double lol2 = GQCP::calculateExpectationValue(S, ao_one_rdm);
    double lol3 = GQCP::calculateExpectationValue(h, ao_one_rdm);
    double lol4 = GQCP::calculateExpectationValue(h2, one_rdm);

    double lol5 = GQCP::calculateExpectationValue(g, ao_two_rdm);
    double lol6 = GQCP::calculateExpectationValue(g2, two_rdm);

    double lol7 = 0;
    for (int i = 0; i<K; i++) {
        for (int j = 0; j<K; j++) {
            lol += ao_one_rdm(i,j)*S(j,i);
        }
    }

    for (int i = 0; i<K; i++) {
        for (int j = 0; j<K; j++) {
            for (int k = 0; k<K; k++) {
                for (int l = 0; l<K; l++) {
                    lol7 += ao_two_rdm(i,j,k,l) * g(i,j,k,l)/2;
                }
            }

        }
    }

    std::cout<<std::endl<<trace_energy<<std::endl<<fci_energy;

    double self1g = 0;
    double self1h = 0;
    double self2h = 0;
    double self2g = 0;
    double inter1g = 0;
    double inter1h = 0;
    for (auto i: gto_list) {
        for (auto j: gto_list) {
            self1h += ao_one_rdm(i,j) * h(i,j);
            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    self1g += ao_two_rdm(i,k,j,l) * g(i,k,j,l)/2;
                }
            }
        }
    }


    for (auto i: gto_list3) {
        for (auto j: gto_list3) {
            self2h += ao_one_rdm(i,j) * h(i,j);
            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    self2g += ao_two_rdm(i,k,j,l) * g(i,k,j,l)/2;
                }
            }
        }
    }


    for (auto i: gto_list3) {
        for (auto j: gto_list) {
            inter1h += 2*ao_one_rdm(i,j) * h(j,i);
            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    inter1g += ao_two_rdm(i,k,j,l) * g(i,k,j,l);
                }
            }
        }
    }

    std::cout<<std::endl;
    std::cout<<self1h<<std::endl;
    std::cout<<self2h<<std::endl;
    std::cout<<self1g<<std::endl;
    std::cout<<self2g<<std::endl;
    std::cout<<inter1g<<std::endl;
    std::cout<<inter1h<<std::endl;
    std::cout<<self1h+self2h+self1g+self2g+inter1g+inter1h<<std::endl;


    std::cout<<"O:"<<self1h+self1g + (inter1h + inter1g)/2<<std::endl;


    /**
     *
     */


    std::cout<<std::endl<<"--------------COMP--------------"<<std::endl;


    double selfgmod = 0;
    double selfgmod2 = 0;

    auto gmod = transform_rdm(g, pa, pa);
    auto gmod2 = transform_rdm(g, pa, pb);


    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {

            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    selfgmod += ao_two_rdm(i,k,j,l) * gmod(i,k,j,l)/2;
                }
            }
        }
    }

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {

            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    selfgmod2 += ao_two_rdm(i,k,j,l) * gmod2(i,k,j,l);
                }
            }
        }
    }


    std::cout<<self1g<<std::endl;
    std::cout<<selfgmod<<std::endl;

    std::cout<<inter1g<<std::endl;
    std::cout<<selfgmod2<<std::endl;


}


BOOST_AUTO_TEST_CASE ( decomposition_H2O_STO_3G ) {

    GQCP::Atom N_1 (7, 0.0, 0.0, 0.0);
    GQCP::Atom N_2 (7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Atom> atoms {N_1, N_2};
    GQCP::Molecule N2 (atoms);

    auto aob = GQCP::AOBasis(N2, "STO-3G");

    auto x = aob.calculateKineticIntegrals();
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<x;
    std::cout<<std::endl;


    GQCP::Atom N_12 (7, 0.0, 0.0, 0.0);
    GQCP::Atom N_22 (7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(3));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Atom> atoms2 {N_12, N_22};
    GQCP::Molecule N22 (atoms2);

    auto aob2 = GQCP::AOBasis(N22, "STO-3G");

    auto x2 = aob2.calculateKineticIntegrals();
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<x2;
    std::cout<<std::endl;





    GQCP::Atom N_13 (7, 0.0, 0.0, 0.0);
    std::vector<GQCP::Atom> atoms3 {N_13};
    GQCP::Molecule N3 (atoms3);

    auto aob3 = GQCP::AOBasis(N3, "STO-3G");

    auto x3 = aob3.calculateKineticIntegrals();
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<x3;
    std::cout<<std::endl;

}
