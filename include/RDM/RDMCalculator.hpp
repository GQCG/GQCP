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
#include "FockSpace/ProductFockSpace.hpp"
#include "FockSpace/SelectedFockSpace.hpp"

#include <memory>


namespace GQCP {


/**
 *  A wrapper around the derived RDMBuilders that provides the functionality of the appropriate derived RDMBuilder for a given Fock space at compile- or runtime.
 */
class RDMCalculator {
private:
    std::shared_ptr<GQCP::BaseRDMBuilder> rdm_builder;

public:
    // CONSTRUCTOR
    /**
     *  Allocate a DOCIRDMBuilder
     *
     *  @param fock_space       the DOCI Fock space
     */
    RDMCalculator(const FockSpace& fock_space);

    /**
     *  Allocate a FCIRDMBuilder
     *
     *  @param fock_space       the FCI Fock space
     */
    RDMCalculator(const ProductFockSpace& fock_space);

    /**
     *  Allocate a SelectedRDMBuilder
     *
     *  @param fock_space       the 'selected' Fock space
     */
    RDMCalculator(const SelectedFockSpace& fock_space);

    /**
     *  A run-time constructor allocating the appropriate derived RDMBuilder
     *
     *  @param fock_space       the Fock space on which the RDMBuilder should be based
     */
    RDMCalculator(const BaseFockSpace& fock_space);


    // PUBLIC METHODS
    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 1-RDMs given a coefficient vector
     */
    OneRDMs calculate1RDMs(const Eigen::VectorXd& x);

    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 2-RDMs given a coefficient vector
     */
    TwoRDMs calculate2RDMs(const Eigen::VectorXd& x);
};


}  // namespace GQCP


#endif  // GQCP_RDMCALCULATOR_HPP
