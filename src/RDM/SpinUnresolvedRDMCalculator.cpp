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
#include "RDM/SpinUnresolvedRDMCalculator.hpp"

#include "RDM/DOCIRDMBuilder.hpp"
#include "RDM/FCIRDMBuilder.hpp"
#include "RDM/SelectedRDMBuilder.hpp"
#include "RDM/SpinUnresolvedFCIRDMBuilder.hpp"



namespace GQCP {


/*
 *  CONSTRUCTOR
 */

/**
 *  Allocate a SpinUnresolvedFCIRDMBuilder
 *
 *  @param fock_space       the Fock space
 */
SpinUnresolvedRDMCalculator::SpinUnresolvedRDMCalculator(const FockSpace& fock_space) :
    rdm_builder(SpinUnresolvedFCIRDMBuilder(fock_space))
{}

/**
 *  A run-time constructor allocating the appropriate derived RDMBuilder and coefficient vector
 *
 *  @param wavefunction       the wave function holding the coefficient vector and a Fock space on which the RDMBuilder should be based
 */
SpinUnresolvedRDMCalculator::SpinUnresolvedRDMCalculator(const SpinUnresolvedWaveFunction& wavefunction) :
    SpinUnresolvedRDMCalculator(dynamic_cast<const FockSpace&>(wavefunction.get_fock_space()))
{
    this->set_coefficients(wavefunction.get_coefficients());
}

/*
 *  PUBLIC METHODS
 */

/**
 *  @return the 1-RDM if a given coefficient vector is set
 */
OneRDM SpinUnresolvedRDMCalculator::calculate1RDM() const {
    if (this->coefficients.rows() == 0) { throw std::logic_error("No vector has been set."); }
    return rdm_builder.calculate1RDM(this->coefficients);
}

/**
 *  @return the 2-RDM if a given coefficient vector is set
 */
TwoRDM SpinUnresolvedRDMCalculator::calculate2RDM() const {
    if (this->coefficients.rows() == 0) { throw std::logic_error("No vector has been set."); }
    return rdm_builder.calculate2RDM(this->coefficients);
}


/**
 *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
 *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
 *
 *  @return an element of the N-RDM, as specified by the given bra and ket indices
 *
 *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
 */
double SpinUnresolvedRDMCalculator::calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices) const {
    if (this->coefficients.rows() == 0) { throw std::logic_error("No vector has been set."); }
    return this->rdm_builder.calculateElement(bra_indices, ket_indices, this->coefficients);
}


}  // namespace GQCP
