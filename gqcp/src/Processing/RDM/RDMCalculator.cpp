// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
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
#include "Processing/RDM/RDMCalculator.hpp"

#include "Processing/RDM/DOCIRDMBuilder.hpp"
#include "Processing/RDM/FCIRDMBuilder.hpp"
#include "Processing/RDM/FrozenCoreDOCIRDMBuilder.hpp"
#include "Processing/RDM/FrozenCoreFCIRDMBuilder.hpp"
#include "Processing/RDM/SelectedRDMBuilder.hpp"


namespace GQCP {


/*
 *  CONSTRUCTOR
 */

/**
 *  Allocate a DOCIRDMBuilder
 *
 *  @param fock_space       the DOCI ONV basis
 */
RDMCalculator::RDMCalculator(const SpinUnresolvedONVBasis& fock_space) :
    rdm_builder (std::make_shared<DOCIRDMBuilder>(fock_space))
{}


/**
 *  Allocate a FCIRDMBuilder
 *
 *  @param fock_space       the FCI ONV basis
 */
RDMCalculator::RDMCalculator(const SpinResolvedONVBasis& fock_space) :
    rdm_builder (std::make_shared<FCIRDMBuilder>(fock_space))
{}


/**
 *  Allocate a SelectedRDMBuilder
 *
 *  @param fock_space       the 'selected' ONV basis
 */
RDMCalculator::RDMCalculator(const SpinResolvedSelectedONVBasis& fock_space) :
    rdm_builder (std::make_shared<SelectedRDMBuilder>(fock_space))
{}


/**
 *  A run-time constructor allocating the appropriate derived RDMBuilder
 *
 *  @param fock_space       the ONV basis on which the RDMBuilder should be based
 */
RDMCalculator::RDMCalculator(const BaseONVBasis& fock_space) {

    switch (fock_space.get_type()){

        case ONVBasisType::SpinUnresolvedONVBasis: {
            this->rdm_builder = std::make_shared<DOCIRDMBuilder>(dynamic_cast<const SpinUnresolvedONVBasis&>(fock_space));
            break;
        }

        case ONVBasisType::SpinResolvedONVBasis: {
            this->rdm_builder = std::make_shared<FCIRDMBuilder>(dynamic_cast<const SpinResolvedONVBasis&>(fock_space));

            break;
        }

        case ONVBasisType::SpinResolvedSelectedONVBasis: {
            this->rdm_builder = std::make_shared<SelectedRDMBuilder>(dynamic_cast<const SpinResolvedSelectedONVBasis&>(fock_space));

            break;
        }

        case ONVBasisType::SpinUnresolvedFrozenONVBasis: {
            this->rdm_builder = std::make_shared<FrozenCoreDOCIRDMBuilder>(dynamic_cast<const SpinUnresolvedFrozenONVBasis&>(fock_space));

            break;
        }

        case ONVBasisType::SpinResolvedFrozenONVBasis: {
            this->rdm_builder = std::make_shared<FrozenCoreFCIRDMBuilder>(dynamic_cast<const SpinResolvedFrozenONVBasis&>(fock_space));

            break;
        }
    }
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return all 1-RDMs if a given coefficient vector is set
 */
OneRDMs<double> RDMCalculator::calculate1RDMs() const {

    if (this->coefficients.rows() == 0) {
        throw std::logic_error("RDMCalculator::calculate1RDMs(): No vector has been set.");
    }

    return rdm_builder->calculate1RDMs(this->coefficients);
}


/**
 *  @return all 2-RDMs if a given coefficient vector is set
 */
TwoRDMs<double> RDMCalculator::calculate2RDMs() const {

    if (this->coefficients.rows() == 0) {
        throw std::logic_error("RDMCalculator::calculate2RDMs(): No vector has been set.");
    }

    return rdm_builder->calculate2RDMs(this->coefficients);
}


/**
 *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
 *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
 *
 *  @return an element of the N-RDM, as specified by the given bra and ket indices
 *
 *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
 */
double RDMCalculator::calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices) const {

    if (this->coefficients.rows() == 0) {
        throw std::logic_error("RDMCalculator::calculateElement(std::vector<size_t>, std::vector<size_t>): No vector has been set.");
    }

    return this->rdm_builder->calculateElement(bra_indices, ket_indices, this->coefficients);
}


}  // namespace GQCP
