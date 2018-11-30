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
#ifndef GQCP_SELECTEDFOCKSPACE_HPP
#define GQCP_SELECTEDFOCKSPACE_HPP


#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "Configuration.hpp"


namespace GQCP {


/**
 *  A class that represents a Fock space that is flexible in the number of states that span it
 *
 *  Configurations are represented as a Configuration: a combination of an alpha and a beta ONV
 */
class SelectedFockSpace : public GQCP::BaseFockSpace {
private:
    size_t N_alpha;  // number of alpha electrons
    size_t N_beta;  // number of beta electrons

    std::vector<GQCP::Configuration> configurations;

    /**
     *  @param onv1     the alpha ONV as a string representation read from right to left
     *  @param onv2     the beta ONV as a string representation read from right to left
     *
     *  @return the configuration that holds both ONVs
     *
     *  IMPORTANT: only works for up to 64 bits!
     */
    Configuration makeConfiguration(const std::string& onv1, const std::string& onv2);

public:
    // CONSTRUCTORS
    SelectedFockSpace() = default;

    /**
     *  A constructor with initial Fock space dimension of 0
     *
     *  @param K            the number of orbitals
     *  @param N_alpha      the number of alpha electrons
     *  @param N_beta       the number of beta electrons
     */
    SelectedFockSpace(size_t K, size_t N_alpha, size_t N_beta);

    /**
     *  A constructor that generates the configurations based off the given ProductFockSpace.
     *
     *  @param fock_space       the ProductFockSpace from which the configurations should be generated
     */
    explicit SelectedFockSpace(const ProductFockSpace& fock_space);

    /**
     *  A constructor that generates the configurations based off the given FockSpace.
     *
     *  @param fock_space       the FockSpace from which the configurations should be generated
     */
    explicit SelectedFockSpace(const FockSpace& fock_space);


    // GETTERS
    size_t get_N_alpha() const { return this->N_alpha; }
    size_t get_N_beta() const { return this->N_beta; }
    const Configuration& get_configuration(size_t index) const { return this->configurations[index]; }
    FockSpaceType get_type() const override { return FockSpaceType::SelectedFockSpace; }


    // PUBLIC METHODS
    /**
     *  Make a configuration (see makeConfiguration()) and add it to this Fock space
     *
     *  @param onv1     the alpha ONV as a string representation read from right to left
     *  @param onv2     the beta ONV as a string representation read from right to left
     */
    void addConfiguration(const std::string& onv1, const std::string& onv2);

    /**
     *  Make configurations (see makeConfiguration()) and add them to the Fock space
     *
     *  @param onv1s     the alpha ONVs as string representations read from right to left
     *  @param onv2s     the beta ONVs as string representations read from right to left
     */
    void addConfiguration(const std::vector<std::string>& onv1s, const std::vector<std::string>& onv2s);
};


}  // namespace GQCP


#endif  // GQCP_SELECTEDFOCKSPACE_HPP
