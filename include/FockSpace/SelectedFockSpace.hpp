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
#ifndef GQCP_SELECTEDSPACEPRODUCT_HPP
#define GQCP_SELECTEDSPACEPRODUCT_HPP


#include "FockSpace/BaseFockSpace.hpp"
#include "Configuration.hpp"

#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace GQCP {


/**
 *  A Fock space for a given set of orbitals and number of alpha and beta electrons.
 *  Where the considered configurations are manually selected and represented as an ONV pair :
 *  a combination of two ONVs, one holding the alpha configuration and holding the bea configuration.
 */
class SelectedFockSpace : public GQCG::BaseFockSpace {
private:
    const size_t N_alpha;  // number of alpha electrons
    const size_t N_beta;  // number of beta electrons

    std::vector<Configuration> selection;  // the selected configurations

    /**
     *  Given @param onv1 and onv2, make a configuration
     */
    Configuration makeConfiguration(std::string onv1, std::string onv2);

public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param K (spatial orbitals), N_alpha and N_beta (electrons);
     *  the initial dimension of the space is 0 as no selections are made.
     */
    SelectedFockSpace(size_t K, size_t N_alpha, size_t N_beta);


    // DESTRUCTORS
    ~SelectedFockSpace() override = default;


    // GETTERS
    size_t get_N_alpha() const { return this->N_alpha; }
    size_t get_N_beta() const { return this->N_beta; }
    Configuration get_Configuration(size_t index) const { return this->selection[index]; }
    FockSpaceType get_type() const override { return FockSpaceType::SelectedFockSpace; }


    // PUBLIC METHODS
    /**
     *  Add configuration(s) (@param onv1(s) and onv2(s)) to the current selected Fock space
     *  increasing its selected dimension.
     */
    void addConfiguration(std::string onv1, std::string onv2);
    void addConfiguration(std::vector<std::string> onv1s, std::vector<std::string> onv2s);
};


}  // namespace GQCP


#endif  // GQCP_SELECTEDSPACEPRODUCT_HPP
