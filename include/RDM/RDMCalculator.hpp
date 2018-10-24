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
#ifndef GQCP_RDMCALCULATOR_HPP
#define GQCP_RDMCALCULATOR_HPP

#include "RDM/BaseRDMBuilder.hpp"
#include "FockSpace/FockSpace.hpp"
#include "FockSpace/FockSpaceProduct.hpp"
#include "FockSpace/SelectedFockSpace.hpp"

#include <memory>


namespace GQCP {


/**
 *  RDMCalculator is a wrapper around the derived RDMBuilders that provides the functionality of the correct derived RDMBuilder
 *  for a given Fock space at compile- or runtime.
 */
class RDMCalculator {
private:
    std::shared_ptr<GQCP::BaseRDMBuilder> rdm_builder;

public:
    // CONSTRUCTOR
    /**
     *  Allocates a derived RDMBuilder based on the nature of the given @param fock_space
     */
    RDMCalculator(const FockSpace& fock_space);  // DOCIRDMBuilder
    RDMCalculator(const FockSpaceProduct& fock_space);  // FCIRDMBuilder
    RDMCalculator(const SelectedFockSpace& fock_space);  // SelectedRDMBuilder
    RDMCalculator(const BaseFockSpace& fock_space);  // runtime


    // PUBLIC METHODS
    /**
     *  @return all 1-RDMs from a coefficient vector @param x
     */
    OneRDMs calculate1RDMs(const Eigen::VectorXd& x);

    /**
     *  @return all 2-RDMs from a coefficient vector @param x
     */
    TwoRDMs calculate2RDMs(const Eigen::VectorXd& x);
};


}  // namespace GQCP


#endif  // GQCP_RDMCALCULATOR_HPP
