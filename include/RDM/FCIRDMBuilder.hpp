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
#ifndef GQCP_FCIRDMBUILDER_HPP
#define GQCP_FCIRDMBUILDER_HPP


#include "FockSpace/ProductFockSpace.hpp"
#include "RDM/BaseRDMBuilder.hpp"
#include "RDM/RDMs.hpp"


namespace GQCP {


/**
 *  A class capable of calculating 1- and 2-RDMs from wave functions expanded in the full CI product Fock space
 */
class FCIRDMBuilder : public BaseRDMBuilder {
    ProductFockSpace fock_space;  // Fock space containing the alpha and beta Fock space


public:
    // CONSTRUCTORS
    explicit FCIRDMBuilder(const ProductFockSpace& fock_space);


    // DESTRUCTOR
    ~FCIRDMBuilder() = default;


    // OVERRIDDEN GETTERS
    BaseFockSpace* get_fock_space() override { return &fock_space; }


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param x        the coefficient vector representing the FCI wave function
     *
     *  @return all 1-RDMs given a coefficient vector
     */
    OneRDMs calculate1RDMs(const Eigen::VectorXd& x) const override;

    /**
     *  @param x        the coefficient vector representing the FCI wave function
     *
     *  @return all 2-RDMs given a coefficient vector
     */
    TwoRDMs calculate2RDMs(const Eigen::VectorXd& x) const override;
};


}  // namespace GQCP


#endif  // GQCP_FCIRDMBUILDER_HPP
