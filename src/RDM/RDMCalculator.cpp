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
 *  Allocate a DOCIRDMBuilder
 *
 *  @param fock_space       the DOCI Fock space
 */
RDMCalculator::RDMCalculator(const FockSpace& fock_space) :
    rdm_builder (std::make_shared<DOCIRDMBuilder>(fock_space))
{}


/**
 *  Allocate a FCIRDMBuilder
 *
 *  @param fock_space       the FCI Fock space
 */
RDMCalculator::RDMCalculator(const ProductFockSpace& fock_space) :
    rdm_builder (std::make_shared<FCIRDMBuilder>(fock_space))
{}


/**
 *  Allocate a SelectedRDMBuilder
 *
 *  @param fock_space       the 'selected' Fock space
 */
RDMCalculator::RDMCalculator(const SelectedFockSpace& fock_space) :
    rdm_builder (std::make_shared<SelectedRDMBuilder>(fock_space))
{}


/**
 *  A run-time constructor allocating the appropriate derived RDMBuilder
 *
 *  @param fock_space       the Fock space on which the RDMBuilder should be based
 */
RDMCalculator::RDMCalculator(const BaseFockSpace& fock_space) {

    switch (fock_space.get_type()){

        case FockSpaceType::FockSpace: {
            this->rdm_builder = std::make_shared<DOCIRDMBuilder>(dynamic_cast<const FockSpace&>(fock_space));

            break;
        }

        case FockSpaceType::ProductFockSpace: {
            this->rdm_builder = std::make_shared<FCIRDMBuilder>(dynamic_cast<const ProductFockSpace&>(fock_space));

            break;
        }

        case FockSpaceType::SelectedFockSpace: {
            this->rdm_builder = std::make_shared<SelectedRDMBuilder>(dynamic_cast<const SelectedFockSpace&>(fock_space));

            break;
        }
    }

}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param x        the coefficient vector representing the wave function
 *
 *  @return all 1-RDMs given a coefficient vector
 */
OneRDMs RDMCalculator::calculate1RDMs(const Eigen::VectorXd& x) const {
    return rdm_builder->calculate1RDMs(x);
}


/**
 *  @param x        the coefficient vector representing the wave function
 *
 *  @return all 2-RDMs given a coefficient vector
 */
TwoRDMs RDMCalculator::calculate2RDMs(const Eigen::VectorXd& x) const {
    return rdm_builder->calculate2RDMs(x);
}


}  // namespace GQCP
