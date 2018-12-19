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
#include "RDM/UnresolvedCIRDMBuilder.hpp"



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


/**
 *  A run-time constructor allocating the appropriate derived RDMBuilder and coefficient vector
 *
 *  @param wavefunction       the wave function holding the coefficient vector and a Fock space on which the RDMBuilder should be based
 */
RDMCalculator::RDMCalculator(const WaveFunction& wavefunction) :
        RDMCalculator(wavefunction.get_fock_space())
{
    this->set_coefficients(wavefunction.get_coefficients());
}

/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  A run-time constructor allocating the appropriate derived Unresolved RDMBuilder and coefficient vector
 *
 *  @param wavefunction       the wave function holding the coefficient vector and a Fock space on which the RDMBuilder should be based
 */
RDMCalculator RDMCalculator::SpinUnresolved(const WaveFunction& wavefunction) {

    RDMCalculator rdm_calculator = RDMCalculator::SpinUnresolved(wavefunction.get_fock_space());
    rdm_calculator.set_coefficients(wavefunction.get_coefficients());
    return rdm_calculator;
};


/**
 *  Allocate an UnresolvedCIRDMBuilder
 *
 *  @param fock_space       the Fock space
 */
RDMCalculator RDMCalculator::SpinUnresolved(const FockSpace& fock_space) {
    RDMCalculator rdm_calculator;
    rdm_calculator.rdm_builder = std::make_shared<UnresolvedCIRDMBuilder>(fock_space);
    return rdm_calculator;
}


/**
 *  A run-time constructor allocating the appropriate derived Unresolved RDMBuilder
 *
 *  @param fock_space       the Fock space on which the RDMBuilder should be based
 */
RDMCalculator RDMCalculator::SpinUnresolved(const BaseFockSpace& fock_space)
{
    switch (fock_space.get_type()){

        case FockSpaceType::FockSpace: {
            return RDMCalculator::SpinUnresolved(dynamic_cast<const FockSpace&>(fock_space));
        }

        case FockSpaceType::ProductFockSpace: {
            throw std::runtime_error ("No ProductFockSpace based implementation for Spin Unresolved RDMS");
        }

        case FockSpaceType::SelectedFockSpace: {
            throw std::runtime_error ("No ProductFockSpace based implementation for Spin Unresolved RDMS");
        }
    }
}

/*
 *  PUBLIC METHODS
 */

/**
 *  @return all 1-RDMs if a given coefficient vector is set
 */
OneRDMs RDMCalculator::calculate1RDMs() const {
    if (this->coefficients.rows() == 0) { throw std::logic_error("No vector has been set."); }
    return rdm_builder->calculate1RDMs(this->coefficients);
}


/**
 *  @return all 2-RDMs if a given coefficient vector is set
 */
TwoRDMs RDMCalculator::calculate2RDMs() const {
    if (this->coefficients.rows() == 0) { throw std::logic_error("No vector has been set."); }
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
    return this->rdm_builder->calculateElement(bra_indices, ket_indices, this->coefficients);
}


}  // namespace GQCP
