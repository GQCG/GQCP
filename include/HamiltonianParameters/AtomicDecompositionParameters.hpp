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
#ifndef GQCP_ATOMICDECOMPOSITIONPARAMETERS_HPP
#define GQCP_ATOMICDECOMPOSITIONPARAMETERS_HPP


#include <Molecule.hpp>
#include "HamiltonianParameters/HamiltonianParameters.hpp"

namespace GQCP {


/**
 *  A struct that holds the atomic decomposed Hamiltonian Parameters.
 *
 *  @tparam Scalar      the scalar type
 */
struct AtomicDecompositionParameters {

    // Molecule molecule;  // decomposed molecule

    std::vector<HamiltonianParameters<double>> net_atomic_parameters;  // vector of net atomic Hamiltonian parameters, E_AA, E_BB
    std::vector<HamiltonianParameters<double>> interaction_parameters;  // vector of interaction Hamiltonian parameters, E_AB
    std::vector<HamiltonianParameters<double>> fragment_parameters;  // vector of total atomic or fragment contributions, E_AA + E_AB/2, E_BB + E_AB/2





    AtomicDecompositionParameters (const HamiltonianParameters<double>& mol_ham_par){

        if (!mol_ham_par.get_ao_basis()) {
            throw std::invalid_argument("HamiltonianParameters::atomicDecomposition(): The Hamiltonian parameters have no underlying AO basis, atomic decomposition is not possible.");
        }

        HamiltonianParameters<double> ham_par = mol_ham_par;
        ham_par.transform<double>(mol_ham_par.get_T_total().inverse());

        auto atoms = mol_ham_par.get_ao_basis()->get_shell_set().atoms();

        if (atoms.size() > 2) {
            throw std::invalid_argument("HamiltonianParameters::atomicDecomposition(): The Hamiltonian parameters are set up for more than 2 atoms, currently only available for diatomic molecules");
        }

        Molecule atom_a ({atoms[0]});
        Molecule atom_b ({atoms[1]});

        auto basisset_name = mol_ham_par.get_ao_basis()->get_shell_set().get_basisset_name();
        AOBasis aobasis_a (atom_a, mol_ham_par.get_ao_basis()->get_shell_set().get_basisset_name());
        AOBasis aobasis_b (atom_b, mol_ham_par.get_ao_basis()->get_shell_set().get_basisset_name());

        auto Ka = aobasis_a.numberOfBasisFunctions();
        auto Kb = aobasis_b.numberOfBasisFunctions();
        auto K = mol_ham_par.get_K();

        auto p_a = SquareMatrix<double>::PartitionMatrix(0, Ka, mol_ham_par.get_K());
        auto p_b = SquareMatrix<double>::PartitionMatrix(Ka, Kb, mol_ham_par.get_K());


        OneElectronOperator<double> h_a = p_a * ham_par.get_h() * p_a;
        OneElectronOperator<double> h_b = p_b * ham_par.get_h() * p_b;
        OneElectronOperator<double> h_ab = p_b * ham_par.get_h() * p_a + p_a * ham_par.get_h() * p_b;

        auto g_a = ham_par.get_g();
        auto g_b = ham_par.get_g();
        auto g_ab = ham_par.get_g();
        auto g_ba = ham_par.get_g();

        g_a.fourModeMultiplication<double>(p_a, SquareMatrix<double>::Identity(K, K), p_a, SquareMatrix<double>::Identity(K, K));
        g_b.fourModeMultiplication<double>(p_b, SquareMatrix<double>::Identity(K, K), p_b, SquareMatrix<double>::Identity(K, K));
        g_ab.fourModeMultiplication<double>(p_a, SquareMatrix<double>::Identity(K, K), p_b, SquareMatrix<double>::Identity(K, K));
        g_ba.fourModeMultiplication<double>(p_b, SquareMatrix<double>::Identity(K, K), p_a, SquareMatrix<double>::Identity(K, K));
        GQCP::TwoElectronOperator<double> g_abba = g_ab.Eigen() + g_ba.Eigen();

        HamiltonianParameters<double> HAA(mol_ham_par.get_ao_basis(), ham_par.get_S(), h_a, g_a, SquareMatrix<double>::Identity(mol_ham_par.get_K(), mol_ham_par.get_K()));
        HamiltonianParameters<double> HBB(mol_ham_par.get_ao_basis(), ham_par.get_S(), h_b, g_b, SquareMatrix<double>::Identity(mol_ham_par.get_K(), mol_ham_par.get_K()));
        HamiltonianParameters<double> HAB(mol_ham_par.get_ao_basis(), ham_par.get_S(), h_ab, g_abba, SquareMatrix<double>::Identity(mol_ham_par.get_K(), mol_ham_par.get_K()), mol_ham_par.get_scalar());
        HamiltonianParameters<double> HA(mol_ham_par.get_ao_basis(), ham_par.get_S(), h_a + h_ab/2, g_a.Eigen() + (0.5)*g_abba.Eigen(), SquareMatrix<double>::Identity(mol_ham_par.get_K(), mol_ham_par.get_K()));
        HamiltonianParameters<double> HB(mol_ham_par.get_ao_basis(), ham_par.get_S(), h_b + h_ab/2, g_b.Eigen() + (0.5)*g_abba.Eigen(), SquareMatrix<double>::Identity(mol_ham_par.get_K(), mol_ham_par.get_K()));

        this->net_atomic_parameters = {HAA, HBB};
        this->interaction_parameters = {HAB};
        this->fragment_parameters = {HA, HB};
    };

    static AtomicDecompositionParameters ALT(const HamiltonianParameters<double>& mol_ham_par) {
            if (!this->get_ao_basis()) {
                    throw std::invalid_argument("HamiltonianParameters::atomicDecomposition(): The Hamiltonian parameters have no underlying AO basis, atomic decomposition is not possible.");
            }

            auto ham_par = *this;
            ham_par.transform(this->T_total.inverse());

            auto atoms = this->get_ao_basis().get_shell_set().atoms();

            if (atoms.size() > 2) {
                    throw std::invalid_argument("HamiltonianParameters::atomicDecomposition(): The Hamiltonian parameters are set up for more than 2 atoms, currently only available for diatomic molecules");
            }

            Molecule atom_a ({atoms[0]});
            Molecule atom_b ({atoms[1]});

            auto basisset_name = this->get_ao_basis().get_shell_set().get_basisset_name();
            AOBasis aobasis_a (atom_a, this->get_ao_basis().get_shell_set().get_basisset_name());
            AOBasis aobasis_b (atom_b, this->get_ao_basis().get_shell_set().get_basisset_name());

            auto Ka = aobasis_a.numberOfBasisFunctions();
            auto Kb = aobasis_b.numberOfBasisFunctions();

            auto p_a = SquareMatrix::PartitionMatrix(0, Ka, this->K);
            auto p_b = SquareMatrix::PartitionMatrix(Ka, Kb, this->K);

            auto v_a = aobasis_a.calculateNuclearIntegrals();
            auto v_b = aobasis_b.calculateNuclearIntegrals();

            auto t_a = aobasis_a.calculateKineticIntegrals();
            auto t_b = aobasis_b.calculateKineticIntegrals();

            auto t_ab = this->ao_basis.calculateKineticIntegrals() - t_a - t_b;
            auto v_ab = this->ao_basis.calculateNuclearIntegrals() - v_a - v_b;

            OneElectronOperator<double> h_a = OneElectronOperator<double>::Zero(this->K, this->K);
            OneElectronOperator<double> h_b = OneElectronOperator<double>::Zero(this->K, this->K);

            h_a.block(0,0,Ka,Ka) = t_a + v_a;
            h_b.block(Ka,Ka,Kb,Kb) = t_b + v_b;

            auto g_a = ham_par.get_g();
            auto g_b = ham_par.get_g();
            auto g_ab = ham_par.get_g();

            g_a.fourModeMultiplication(p_a, SquareMatrix::Identiy(this->K, this->K), p_a, SquareMatrix::Identiy(this->K, this->K));
            g_b.fourModeMultiplication(p_b, SquareMatrix::Identiy(this->K, this->K), p_b, SquareMatrix::Identiy(this->K, this->K));
            g_ab.fourModeMultiplication(p_a, SquareMatrix::Identiy(this->K, this->K), p_b, SquareMatrix::Identiy(this->K, this->K));

            HamiltonianParameters<double> HAA(this->ao_basis, ham_par.get_S(), h_a, g_a, SquareMatrix::Identiy(this->K, this->K));
            HamiltonianParameters<double> HBB(this->ao_basis, ham_par.get_S(), h_b, g_b, SquareMatrix::Identiy(this->K, this->K));
            HamiltonianParameters<double> HAB(this->ao_basis, ham_par.get_S(), ham_par.get_h()-h_a-h_b, 2*g_ab, SquareMatrix::Identiy(this->K, this->K));
            HamiltonianParameters<double> HA(this->ao_basis, ham_par.get_S(), (ham_par.get_h()-h_a-h_b)/2, g_ab, SquareMatrix::Identiy(this->K, this->K));
            //HamiltonianParameters<double> HA();
            //HamiltonianParameters<double> HB;
    };
};

}  // namespace GQCP


#endif //GQCP_ATOMICDECOMPOSITIONPARAMETERS_HPP
