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
 *  Given @param onv1 and onv2, make a configuration
 */
Configuration SelectedFockSpace::makeConfiguration(std::string onv1, std::string onv2){

    boost::dynamic_bitset<> alpha_transfer (onv1);
    boost::dynamic_bitset<> beta_transfer (onv2);

    if (alpha_transfer.size() != this->K | beta_transfer.size() != this->K) {
        throw std::invalid_argument("Added configuration is not compatible with the orbitals number of the Fock space");
    }

    if (alpha_transfer.count() != this->N_alpha | beta_transfer.count() != this->N_beta) {
        throw std::invalid_argument("Added configuration is not compatible with the electron numbers of the Fock space");
    }

    size_t alpha_s = alpha_transfer.to_ulong();
    size_t beta_s = beta_transfer.to_ulong();

    ONV alpha (alpha_transfer.size(), alpha_s);
    ONV beta (beta_transfer.size(), beta_s);

    return Configuration { alpha, beta };
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



/*
 *  PUBLIC METHODS
 */

/**
 *  Add configuration(s) (@param onv1(s) and onv2(s)) to the current selected Fock space
 *  increasing its selected dimension.
 */
void SelectedFockSpace::addConfiguration(std::string onv1, std::string onv2){

    this->dim++;

    Configuration configuration = makeConfiguration(onv1, onv2);
    selection.push_back(configuration);
}

void SelectedFockSpace::addConfiguration(std::vector<std::string> onv1s, std::vector<std::string> onv2s){

    if (onv1s.size() != onv2s.size()) {
        throw std::invalid_argument("Size of both ONV entry vectors do not match");
    }

    for (int i = 0; i<onv1s.size(); i++) {
        this->addConfiguration(onv1s[i], onv2s[i]);
    }
}


}  // namespace GQCP
