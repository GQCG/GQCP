// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <math.h>       // sqrt

#include "LibintCommunicator.hpp"


namespace GQCP {


/**
 *  @param ao_basis     the AO basis in which the molecular Hamiltonian should be expressed
 *
 *  @return Hamiltonian parameters corresponding to the molecular Hamiltonian. The molecular Hamiltonian has
 *      - one-electron contributions:
 *          - kinetic
 *          - nuclear attraction
 *      - two-electron contributions:
 *          - Coulomb repulsion
 */
GQCP::HamiltonianParameters constructMolecularHamiltonianParameters(std::shared_ptr<GQCP::AOBasis> ao_basis) {

    // Calculate the integrals for the molecular Hamiltonian
    auto S = GQCP::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::overlap, *ao_basis);
    auto T = GQCP::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::kinetic, *ao_basis);
    auto V = GQCP::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::nuclear, *ao_basis);
    auto H = T + V;
    
    auto g = GQCP::LibintCommunicator::get().calculateTwoElectronIntegrals(libint2::Operator::coulomb, *ao_basis);
    
    
    // Construct the initial transformation matrix: the identity matrix
    auto nbf = ao_basis->get_number_of_basis_functions();
    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(nbf, nbf);
    
    
    return HamiltonianParameters(ao_basis, S, H, g, C);
}


/**
 *  @param K        the number of orbitals
 *
 *  @return a set of random Hamiltonian parameters with values uniformly distributed between [-1,1]
 */
GQCP::HamiltonianParameters constructRandomHamiltonianParameters(size_t K) {

    GQCP::OneElectronOperator S (Eigen::MatrixXd::Identity(K, K));  // the underlying orbital basis can be chosen as orthonormal, since the form of the underlying orbitals doesn't really matter
    Eigen::MatrixXd C (Eigen::MatrixXd::Identity(K, K));  // the transformation matrix C here doesn't really mean anything, because it doesn't link to any AO basis

    GQCP::OneElectronOperator H (Eigen::MatrixXd::Random(K, K));  // uniformly distributed between [-1,1]


    // Unfortunately, the Tensor module provides uniform random distributions between [0, 1]
    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();

    // Move the distribution from [0, 1] -> [-1, 1]
    for (size_t i = 0; i < K; i++) {
        for (size_t j = 0; j < K; j++) {
            for (size_t k = 0; k < K; k++) {
                for (size_t l = 0; l < K; l++) {
                    g_tensor(i,j,k,l) = 2*g_tensor(i,j,k,l) - 1;  // scale from [0, 1] -> [0, 2] -> [-1, 1]
                }
            }
        }
    }
    GQCP::TwoElectronOperator g (g_tensor);

    std::shared_ptr<GQCP::AOBasis> ao_basis;  // nullptr because it doesn't make sense to set an AOBasis

    return GQCP::HamiltonianParameters(ao_basis, S, H, g, C);
}
    
    
/**
 *  @param fcidump_file     the name of the FCIDUMP file
 *
 *  @return Hamiltonian parameters corresponding to the contents of an FCIDUMP file
 */
GQCP::HamiltonianParameters readFCIDUMPFile(const std::string& fcidump_filename) {
    
    // Find the extension of the given path (https://stackoverflow.com/a/51992)
    std::string extension;
    std::string::size_type idx = fcidump_filename.rfind('.');
    
    if (idx != std::string::npos) {
        extension = fcidump_filename.substr(idx+1);
    } else {
        throw std::runtime_error("I did not find an extension in your given path.");
    }
    
    if (!(extension == "FCIDUMP")) {
        throw std::runtime_error("You did not provide a .FCIDUMP file name");
    }
    
    // If the xyz_filename isn't properly converted into an input file stream, we assume the user supplied a wrong file
    std::ifstream input_file_stream (fcidump_filename);
    
    if (!input_file_stream.good()) {
        throw std::runtime_error("The provided FCIDUMP file is illegible. Maybe you specified a wrong path?");
    }
    
    
    
    // Do the actual parsing
    
    //  Get the number of orbitals to check if it's a valid FCIDUMP file
    std::string start_line;  // first line contains orbitals and electron count
    std::getline(input_file_stream, start_line);
    std::stringstream linestream (start_line);
    
    size_t K = 0;
    char iter;
    
    while (linestream >> iter) {
        if (iter == '=') {
            linestream >> K;  // right here we have the number of orbitals
            break;  // we can finish reading the linestream after we found K
        }
    }
    
    if (K == 0) {
        throw std::invalid_argument("The .FCIDUMP-file is invalid: could not read a number of orbitals.");
    }
    
    
    Eigen::MatrixXd h_SO = Eigen::MatrixXd::Zero(K, K);
    Eigen::Tensor<double, 4> g_SO (K, K, K, K);
    g_SO.setZero();
    
    //  Skip 3 lines
    for (size_t counter = 0; counter < 3; counter++) {
        std::getline(input_file_stream, start_line);
    }
    
    
    //  Start reading in the one- and two-electron integrals
    double x;
    size_t i, j, a, b;
    
    std::string line;
    while (std::getline(input_file_stream, line)) {
        std::istringstream iss (line);
        
        // Based on what the values of the indices are, we can read one-electron integrals, two-electron integrals and the internuclear repulsion energy
        //  See also (http://hande.readthedocs.io/en/latest/manual/integrals.html)
        //  I think the documentation is a bit unclear for the two-electron integrals, but we can rest assured that FCIDUMP files give the two-electron integrals in CHEMIST's notation.
        iss >> x >> i >> a >> j >> b;
        
        //  Internuclear repulsion energy (skipped)
        if ((i == 0) && (j == 0) && (a == 0) && (b == 0)) {}
        
        //  Single-particle eigenvalues (skipped)
        else if ((a == 0) && (j == 0) && (b == 0)) {}
        
        //  One-electron integrals (h_core)
        else if ((j == 0) && (b == 0)) {
            size_t p = i - 1;
            size_t q = a - 1;
            h_SO(p,q) = x;
            
            // Apply the permutational symmetry for real orbitals
            h_SO(q,p) = x;
        }
        
        //  Two-electron integrals are given in CHEMIST'S NOTATION, so just copy them over
        else if ((i > 0) && (a > 0) && (j > 0) && (b > 0)) {
            size_t p = i - 1;
            size_t q = a - 1;
            size_t r = j - 1;
            size_t s = b - 1;
            g_SO(p,q,r,s) = x;
            
            // Apply the permutational symmetries for real orbitals
            g_SO(p,q,s,r) = x;
            g_SO(q,p,r,s) = x;
            g_SO(q,p,s,r) = x;
            
            g_SO(r,s,p,q) = x;
            g_SO(s,r,p,q) = x;
            g_SO(r,s,q,p) = x;
            g_SO(s,r,q,p) = x;
        }
    }  // while loop
    
    
    // Make the ingredients to construct HamiltonianParameters
    std::shared_ptr<GQCP::AOBasis> ao_basis;  // nullptr
    GQCP::OneElectronOperator S (Eigen::MatrixXd::Identity(K, K));
    GQCP::OneElectronOperator H_core (h_SO);
    GQCP::TwoElectronOperator G (g_SO);
    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(K, K);

    return GQCP::HamiltonianParameters(ao_basis, S, H_core, G, C);
}


/**
 *  @param upper_triagonal      an upper triagonal representation of the Hubbard hopping matrix
 *
 *  @return Hubbard Hamiltonian parameters generated from the Hubbard hopping matrix
 */
GQCP::HamiltonianParameters constructHubbardParameters(const Eigen::VectorXd& upper_triagonal) {

    // The dimension of matrix K is related to the length of the triagonal vector :
    // K (K + 1) = 2*X
    // K = (sqrt(1+4X) - 1)/2

    size_t X = upper_triagonal.rows();
    size_t K = (static_cast<size_t>(sqrt(1 + 8*X) - 1))/2;



    if (K * (K+1) != 2*X) {
        throw std::invalid_argument("Passed vector was not the triagonal of a square matrix");
    }

    Eigen::MatrixXd h_SO = Eigen::MatrixXd::Zero(K, K);
    Eigen::Tensor<double, 4> g_SO (K, K, K, K);
    g_SO.setZero();

    size_t triagonal_index = 0;
    for (size_t i = 0; i < K; i++) {
        for (size_t j = i; j < K; j++) {

            if (i == j){
                g_SO(i,i,i,i) = upper_triagonal(triagonal_index);
            } else {
                h_SO (i,j) = upper_triagonal(triagonal_index);
                h_SO (j,i) = upper_triagonal(triagonal_index);
            }
            triagonal_index++;
        }
    }


    // Make the ingredients to construct HamiltonianParameters
    std::shared_ptr<GQCP::AOBasis> ao_basis;  // nullptr
    GQCP::OneElectronOperator S (Eigen::MatrixXd::Identity(K, K));
    GQCP::OneElectronOperator H_core (h_SO);
    GQCP::TwoElectronOperator G (g_SO);
    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(K, K);

    return GQCP::HamiltonianParameters(ao_basis, S, H_core, G, C);
}


}  // namespace GQCP
