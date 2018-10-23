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
#include "FockSpace/FockSpaceProduct.hpp"
#include "Configuration.hpp"

#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace GQCP {


/**
 *  A Fock space that is flexible in the number of states that span it.
 *  The configurations are represented as a Configuration:
 *  a combination of two ONVs, holding the alpha and beta ONVs.
 */
class SelectedFockSpace : public GQCP::BaseFockSpace {
private:
    size_t N_alpha;  // number of alpha electrons
    size_t N_beta;  // number of beta electrons

    std::vector<GQCP::Configuration> configurations;

    /**
     *  Member taking two string arguments and creating a Configuration
     *  @param onv1 a string representation read from right to left
     *  @param onv2 a string representation read from right to left
     *  @return a Configuration
     *  !!! only works for up to 64 bits !!!
     */
    Configuration makeConfiguration(const std::string& onv1, const std::string& onv2);

public:
    // CONSTRUCTORS
    SelectedFockSpace() = default;
    /**
     *  Constructor given a @param K (spatial orbitals), N_alpha and N_beta (electrons);
     *  the initial dimension of the space is 0, as no selections are made.
     */
    SelectedFockSpace(size_t K, size_t N_alpha, size_t N_beta);

    /**
     * Constructor that generates expansion of a given FockSpaceProduct
     * @param fock_space generated Fock space
     */
    explicit SelectedFockSpace(const FockSpaceProduct& fock_space);

    /**
     * Constructor that generates expansion of a given FockSpace
     * @param fock_space generated Fock space
     */
    explicit SelectedFockSpace(const FockSpace& fock_space);


    // GETTERS
    size_t get_N_alpha() const { return this->N_alpha; }
    size_t get_N_beta() const { return this->N_beta; }
    Configuration get_configuration(size_t index) const { return this->configurations[index]; }
    FockSpaceType get_type() const override { return FockSpaceType::SelectedFockSpace; }


    // PUBLIC METHODS
    /**
     *  Member taking two string arguments to add a Configuration
     *  @see makeConfiguration()
     *  add a Configuration to @var configurations
     */
    void addConfiguration(const std::string& onv1, const std::string& onv2);

    /**
     *  Member taking two vector<string> arguments to add Configurations
     *  @see makeConfiguration()
     *  add multiple Configurations to @var configurations
     */
    void addConfiguration(const std::vector<std::string>& onv1s, const std::vector<std::string>& onv2s);
};


}  // namespace GQCP


#endif  // GQCP_SELECTEDFOCKSPACE_HPP
