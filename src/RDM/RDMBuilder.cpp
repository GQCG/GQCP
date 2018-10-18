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
#include "RDM/RDMBuilder.hpp"

#include "RDM/DOCIRDMBuilder.hpp"
#include "RDM/FCIRDMBuilder.hpp"


namespace GQCP {

/*
 *  CONSTRUCTOR
 */

/**
 *  Allocates a DOCIRDMBuilder based on @param fock_space
 */
RDMBuilder::RDMBuilder(const FockSpace& fock_space) {
    rdm_builder = std::make_shared<GQCG::DOCIRDMBuilder>(fock_space);
}


/**
 *  Allocates a FCIRDMBuilder based on @param fock_space
 */
RDMBuilder::RDMBuilder(const FockSpaceProduct& fock_space) {
    rdm_builder = std::make_shared<GQCG::FCIRDMBuilder>(fock_space);
}


/**
 *  Allocates the correct derived BaseRDMBuilder based on @param fock_space
 */
RDMBuilder::RDMBuilder(const BaseFockSpace& fock_space) {

    FockSpaceType fock_space_type = fock_space.get_fock_space_type();

    switch (fock_space_type){

        case FockSpaceType::FockSpace: {
            rdm_builder = std::make_shared<GQCG::DOCIRDMBuilder>(dynamic_cast<const GQCG::FockSpace&>(fock_space));
        }

        case FockSpaceType::FockSpaceProduct: {
            rdm_builder = std::make_shared<GQCG::FCIRDMBuilder>(dynamic_cast<const GQCG::FockSpaceProduct&>(fock_space));
        }
    }

}


}  // namespace GQCP