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
#include "FockSpace/SelectedFockSpace.hpp"


#include "boost/dynamic_bitset.hpp"


namespace GQCP {


/*
 *  PRIVATE METHODS
 */

/**
 *  Member taking two string arguments and creating a Configuration
 *  @param onv1 a string representation read from right to left
 *  @param onv2 a string representation read from right to left
 *  @return a Configuration
 *  !!! only works for up to 64 bits !!!
 */
Configuration SelectedFockSpace::makeConfiguration(const std::string& onv1, const std::string& onv2){

    boost::dynamic_bitset<> alpha_transfer (onv1);
    boost::dynamic_bitset<> beta_transfer (onv2);

    if (alpha_transfer.size() != this->K | beta_transfer.size() != this->K) {
        throw std::invalid_argument("Given string representations for ONVs are not compatible with the number of orbitals of the Fock space");
    }

    if (alpha_transfer.count() != this->N_alpha | beta_transfer.count() != this->N_beta) {
        throw std::invalid_argument("Given string representations for ONVs are not compatible with the number of orbitals of the Fock space");
    }

    size_t alpha_s = alpha_transfer.to_ulong();
    size_t beta_s = beta_transfer.to_ulong();

    ONV alpha (alpha_transfer.size(), alpha_s);
    ONV beta (beta_transfer.size(), beta_s);

    return Configuration {alpha, beta};
}



/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor given a @param K (spatial orbitals), N_alpha and N_beta (electrons);
 *  the initial dimension of the space is 0 as no selections are made.
 */
SelectedFockSpace::SelectedFockSpace(size_t K, size_t N_alpha, size_t N_beta) :
        BaseFockSpace(K, 0),
        N_alpha (N_alpha),
        N_beta (N_beta)
{}

/**
 * Constructor that generates expansion of a given FockSpaceProduct
 * @param fock_space generated Fock space
 */
SelectedFockSpace::SelectedFockSpace(const FockSpaceProduct& fock_space) :
    SelectedFockSpace (fock_space.get_K(), fock_space.get_N_alpha(), fock_space.get_N_beta())
{

    std::vector<Configuration> configurations;

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    ONV alpha = fock_space_alpha.get_ONV(0);
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {

        ONV beta = fock_space_beta.get_ONV(0);
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {

            configurations.push_back(Configuration {alpha, beta});

            if (I_beta < dim_beta - 1) {  // prevent the last permutation to occur
                fock_space_beta.setNext(beta);
            }
        }
        if (I_alpha < dim_alpha - 1) {  // prevent the last permutation to occur
            fock_space_alpha.setNext(alpha);
        }
    }
    this->dim = fock_space.get_dimension();
    this->configurations = configurations;

}

/**
 * Constructor that generates expansion of a given FockSpace
 * @param fock_space generated Fock space
 */
SelectedFockSpace::SelectedFockSpace(const FockSpace& fock_space)  :
        SelectedFockSpace (fock_space.get_K(), fock_space.get_N(), fock_space.get_N())
{

    std::vector<Configuration> configurations;

    // Current workaround to call non-const functions
    FockSpace fock_space_single = fock_space;


    auto dim = fock_space_single.get_dimension();

    // Iterate over the Fock space and add all onvs as doubly occupied configurations
    ONV onv = fock_space_single.get_ONV(0);
    for (size_t I = 0; I < dim; I++) {

        configurations.push_back(Configuration {onv, onv});

        if (I < dim - 1) {  // prevent the last permutation to occur
            fock_space_single.setNext(onv);
        }

    }

    this->dim = dim;
    this->configurations = configurations;

}


/*
 *  PUBLIC METHODS
 */

/**
 *  Member taking two string arguments to add a Configuration
 *  @see makeConfiguration()
 *  add a Configuration to @var configurations
 */
void SelectedFockSpace::addConfiguration(const std::string& onv1, const std::string& onv2){

    this->dim++;

    Configuration configuration = makeConfiguration(onv1, onv2);
    configurations.push_back(configuration);
}

/**
 *  Member taking two vector<string> arguments to add Configurations
 *  @see makeConfiguration()
 *  add multiple Configurations to @var configurations
 */
void SelectedFockSpace::addConfiguration(const std::vector<std::string>& onv1s, const std::vector<std::string>& onv2s){

    if (onv1s.size() != onv2s.size()) {
        throw std::invalid_argument("Size of both ONV entry vectors do not match");
    }

    for (int i = 0; i<onv1s.size(); i++) {
        this->addConfiguration(onv1s[i], onv2s[i]);
    }
}


}  // namespace GQCP
