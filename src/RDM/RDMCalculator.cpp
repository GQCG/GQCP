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
#include "RDM/RDMCalculator.hpp"

#include "RDM/DOCIRDMBuilder.hpp"
#include "RDM/FCIRDMBuilder.hpp"
#include "RDM/SelectedRDMBuilder.hpp"



namespace GQCP {


/*
 *  CONSTRUCTOR
 */

/**
 *  Allocates a DOCIRDMBuilder based on @param fock_space
 */
RDMCalculator::RDMCalculator(const FockSpace& fock_space) {
    rdm_builder = std::make_shared<GQCP::DOCIRDMBuilder>(fock_space);
}


/**
 *  Allocates a FCIRDMBuilder based on @param fock_space
 */
RDMCalculator::RDMCalculator(const FockSpaceProduct& fock_space) {
    rdm_builder = std::make_shared<GQCP::FCIRDMBuilder>(fock_space);
}


/**
 *  Allocates a SelectedRDMBuilder based on @param fock_space
 */
RDMCalculator::RDMCalculator(const SelectedFockSpace& fock_space) {
    rdm_builder = std::make_shared<GQCP::SelectedRDMBuilder>(fock_space);
}


/**
 *  Allocates the correct derived BaseRDMBuilder based on @param fock_space
 */
RDMCalculator::RDMCalculator(const BaseFockSpace& fock_space) {

    switch (fock_space.get_type()){

        case FockSpaceType::FockSpace: {
            rdm_builder = std::make_shared<GQCP::DOCIRDMBuilder>(dynamic_cast<const GQCP::FockSpace&>(fock_space));

            break;
        }

        case FockSpaceType::FockSpaceProduct: {
            rdm_builder = std::make_shared<GQCP::FCIRDMBuilder>(dynamic_cast<const GQCP::FockSpaceProduct&>(fock_space));

            break;
        }

        case FockSpaceType::SelectedFockSpace: {
            rdm_builder = std::make_shared<GQCP::SelectedRDMBuilder>(dynamic_cast<const GQCP::SelectedFockSpace&>(fock_space));

            break;
        }
    }

}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return all 1-RDMs from a coefficient vector @param x
 */
OneRDMs RDMCalculator::calculate1RDMs(const Eigen::VectorXd& x) {
    return rdm_builder->calculate1RDMs(x);
}


/**
 *  @return all 2-RDMs from a coefficient vector @param x
 */
TwoRDMs RDMCalculator::calculate2RDMs(const Eigen::VectorXd& x) {
    return rdm_builder->calculate2RDMs(x);
}


}  // namespace GQCP
