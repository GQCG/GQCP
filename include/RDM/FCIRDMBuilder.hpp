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


#include "FockSpace/FockSpaceProduct.hpp"
#include "RDM/RDMBuilder.hpp"
#include "RDM/RDMs.hpp"


namespace GQCP {


/**
 *  FCIRDMBuilder is a class for the calculation of a density matrix from a given wave function
 *  or coefficient expansion in the Full CI Fock space
 */
class FCIRDMBuilder : public RDMBuilder {
    FockSpaceProduct fock_space;  // fock space containing the alpha and beta Fock space


public:
    // CONSTRUCTOR
    explicit FCIRDMBuilder(const FockSpaceProduct& fock_space);


    // DESTRUCTOR
    ~FCIRDMBuilder() = default;


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return all 1-RDMs from a coefficient vector @param x
     */
    OneRDMs calculate1RDMs(const Eigen::VectorXd& x) override;

    /**
     *  @return all 2-RDMs from a coefficient vector @param x
     */
    TwoRDMs calculate2RDMs(const Eigen::VectorXd& x) override;

    /**
     *  @return the Fock space of the RDMBuilder
     */
    BaseFockSpace* get_fock_space() override { return &fock_space; }
};


}  // namespace GQCP


#endif  // GQCP_FCIRDMBUILDER_HPP
