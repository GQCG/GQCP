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

#include <boost/dynamic_bitset.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace GQCP {


/*
 *  PRIVATE METHODS
 */

/**
 *  @param onv1     the alpha ONV as a string representation read from right to left
 *  @param onv2     the beta ONV as a string representation read from right to left
 *
 *  @return the configuration that holds both ONVs
 *
 *  IMPORTANT: only works for up to 64 bits!
 */
Configuration SelectedFockSpace::makeConfiguration(const std::string& onv1, const std::string& onv2) const {

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

    ONV alpha (this->K, this->N_alpha, alpha_s);
    ONV beta (this->K, this->N_beta, beta_s);

    return Configuration {alpha, beta};
}



/*
 *  CONSTRUCTORS
 */

/**
 *  A constructor with initial Fock space dimension of 0
 *
 *  @param K            the number of orbitals
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 */
SelectedFockSpace::SelectedFockSpace(size_t K, size_t N_alpha, size_t N_beta) :
    BaseFockSpace(K, 0),
    N_alpha (N_alpha),
    N_beta (N_beta)
{}


/**
 *  A constructor that generates the configurations based off the given ProductFockSpace.
 *
 *  @param fock_space       the ProductFockSpace from which the configurations should be generated
 */
SelectedFockSpace::SelectedFockSpace(const ProductFockSpace& fock_space) :
    SelectedFockSpace (fock_space.get_K(), fock_space.get_N_alpha(), fock_space.get_N_beta())
{

    std::vector<Configuration> configurations;

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    ONV alpha = fock_space_alpha.makeONV(0);
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {

        ONV beta = fock_space_beta.makeONV(0);
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {

            configurations.push_back(Configuration {alpha, beta});

            if (I_beta < dim_beta - 1) {  // prevent the last permutation to occur
                fock_space_beta.setNextONV(beta);
            }
        }
        if (I_alpha < dim_alpha - 1) {  // prevent the last permutation to occur
            fock_space_alpha.setNextONV(alpha);
        }
    }
    this->dim = fock_space.get_dimension();
    this->configurations = configurations;

}


/**
 *  A constructor that generates the configurations based off the given FockSpace.
 *
 *  @param fock_space       the FockSpace from which the configurations should be generated
 */
SelectedFockSpace::SelectedFockSpace(const FockSpace& fock_space)  :
        SelectedFockSpace (fock_space.get_K(), fock_space.get_N(), fock_space.get_N())
{

    std::vector<Configuration> configurations;

    auto dim = fock_space.get_dimension();

    // Iterate over the Fock space and add all onvs as doubly occupied configurations
    ONV onv = fock_space.makeONV(0);
    for (size_t I = 0; I < dim; I++) {

        configurations.push_back(Configuration {onv, onv});

        if (I < dim - 1) {  // prevent the last permutation to occur
            fock_space.setNextONV(onv);
        }

    }

    this->dim = dim;
    this->configurations = configurations;

}



/*
 *  ITERATOR
 */

    /*
     *  CONSTRUCTORS
     */

//    /**
//     *  Place an iterator in a possibly invalid state (the ONV does not correspond to the address)
//     *
//     *  @note This constructor is only implemented to provide a proper SelectedFockSpace.end() iterator with an address that is higher than the 'last' ONV in the Fock space
//     *
//     *  @param fock_space       the Fock space that should be iterated over
//     *  @param address          the address of the current ONV
//     *  @param onv              the current ONV
//     */
//    SelectedFockSpace::Iterator::Iterator(const SelectedFockSpace& fock_space, size_t address, const ONV& onv) :
//        fock_space (&fock_space),
//        address (address),
//        onv (onv)
//    {}


    /**
     *  @param fock_space       the Fock space that should be iterated over
     *  @param address          the address of the current Configuration
     */
    SelectedFockSpace::Iterator::Iterator(const SelectedFockSpace& fock_space, size_t address) :
        fock_space (&fock_space),
        address (address)
    {}


    /**
     *  Constructor that starts at the Configuration with address 0
     *
     *  @param fock_space       the Fock space that should be iterated over
     */
    SelectedFockSpace::Iterator::Iterator(const SelectedFockSpace& fock_space) :
        SelectedFockSpace::Iterator(fock_space, 0)
    {}



    /*
     *  OPERATORS
     */

    /**
     *  Move the iterator forward (pre-increment)
     *
     *  @return a reference to the updated iterator
     */
    SelectedFockSpace::Iterator& SelectedFockSpace::Iterator::operator++() {
        this->address++;
        return *this;
    }


    /**
     *  @param other            the other iterator
     *
     *  @return if this iterator is the same as the other
     */
    bool SelectedFockSpace::Iterator::operator==(const Iterator& other) const {
        return this->address == other.address;
    }


    /**
     *  @param other            the other iterator
     *
     *  @return if this iterator is not the same as the other
     */
    bool SelectedFockSpace::Iterator::operator!=(const Iterator& other) const {
        return !(this->operator==(other));
    }


    /**
     *  @return the current Configuration
     */
    const Configuration& SelectedFockSpace::Iterator::operator*() const {
        return this->fock_space->get_configuration(this->address);
    }



    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the current Configuration
     */
    const Configuration& SelectedFockSpace::Iterator::currentConfiguration() const {
        return this->operator*();
    }


    /**
     *  @return the alpha-ONV of the current Configuration
     */
    const ONV& SelectedFockSpace::Iterator::currentAlphaONV() const {
        return this->currentConfiguration().onv_alpha;
    }

    /**
     *  @return the beta-ONV of the current Configuration
     */
    const ONV& SelectedFockSpace::Iterator::currentBetaONV() const {
        return this->currentConfiguration().onv_beta;
    }

    /**
     *  @return the address of the current Configuration
     */
    size_t SelectedFockSpace::Iterator::currentAddress() const {
        return this->address;
    }



/*
 *  PUBLIC METHODS
 */

/**
 *  Make a configuration (see makeConfiguration()) and add it to this Fock space
 *
 *  @param onv1     the alpha ONV as a string representation read from right to left
 *  @param onv2     the beta ONV as a string representation read from right to left
 */
void SelectedFockSpace::addConfiguration(const std::string& onv1, const std::string& onv2) {

    this->dim++;

    Configuration configuration = makeConfiguration(onv1, onv2);
    configurations.push_back(configuration);
}


/**
 *  Make configurations (see makeConfiguration()) and add them to the Fock space
 *
 *  @param onv1s     the alpha ONVs as string representations read from right to left
 *  @param onv2s     the beta ONVs as string representations read from right to left
 */
void SelectedFockSpace::addConfiguration(const std::vector<std::string>& onv1s, const std::vector<std::string>& onv2s){

    if (onv1s.size() != onv2s.size()) {
        throw std::invalid_argument("Size of both ONV entry vectors do not match");
    }

    for (size_t i = 0; i < onv1s.size(); i++) {
        this->addConfiguration(onv1s[i], onv2s[i]);
    }
}


/**
 *  @return an iterator pointing at the first ONV in the Fock space
 */
SelectedFockSpace::Iterator SelectedFockSpace::begin() {
    return SelectedFockSpace::Iterator(*this);
}


/**
 *  @return an iterator pointing at the last ONV in the Fock space
 */
SelectedFockSpace::Iterator SelectedFockSpace::end() {
    size_t address = this->dim;  // this is an 'invalid' address (i.e. there is no ONV corresponding to this address), but the Iterator design requires .end() to be an iterator following the last element
    return SelectedFockSpace::Iterator(*this, address);
}


}  // namespace GQCP
