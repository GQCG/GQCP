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
#ifndef GQCP_RDMBUILDER_HPP
#define GQCP_RDMBUILDER_HPP

#include "RDM/BaseRDMBuilder.hpp"
#include "FockSpace/FockSpace.hpp"
#include "FockSpace/FockSpaceProduct.hpp"

#include <memory>


namespace GQCP {


/**
 *  RDMBuilder is a wrapper around the RDMBuilders that provides the functionality of the correct RDMBuilder
 *  for a given Fock space at compile- or runtime.
 */


class RDMBuilder {
private:
    std::shared_ptr<GQCG::BaseRDMBuilder> rdm_builder;
public:
    // CONSTRUCTOR
    /**
     *  Allocates a DOCIRDMBuilder based on @param fock_space
     */
    RDMBuilder(const FockSpace& fock_space);

    /**
     *  Allocates a FCIRDMBuilder based on @param fock_space
     */
    RDMBuilder(const FockSpaceProduct& fock_space);

    /**
     *  Allocates the correct derived BaseRDMBuilder based on @param fock_space
     */
    RDMBuilder(const BaseFockSpace& fock_space);


    // PUBLIC METHODS
    /**
     *  @return all 1-RDMs from a coefficient vector @param x
     */
    OneRDMs calculate1RDMs(const Eigen::VectorXd &x) { return rdm_builder->calculate1RDMs(x); }

    /**
     *  @return all 2-RDMs from a coefficient vector @param x
     */
    TwoRDMs calculate2RDMs(const Eigen::VectorXd &x) { return rdm_builder->calculate2RDMs(x); }
};


}  // namespace GQCP


#endif  // GQCP_RDMBUILDER_HPP
