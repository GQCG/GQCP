// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/CoulombRepulsionOperator.hpp"
#include "Operator/FirstQuantized/KineticOperator.hpp"
#include "Operator/FirstQuantized/NuclearAttractionOperator.hpp"


namespace GQCP {


/**
 *  The first-quantized, molecular electronic Hamiltonian.
 */
class FQMolecularHamiltonian {
protected:
    // The kinetic energy operator.
    KineticOperator T;

    // The nuclear attraction operator.
    NuclearAttractionOperator V;

    // The two-electron repulsion operator.
    CoulombRepulsionOperator g;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param T            The kinetic energy operator.
     *  @param V            The nuclear attraction operator.
     *  @param g            The two-electron repulsion operator.
     */
    FQMolecularHamiltonian(const KineticOperator& T, const NuclearAttractionOperator& V, const CoulombRepulsionOperator& g);


    /**
     *  Construct a `FQMolecularHamiltonian` from a `Molecule`.
     * 
     *  @param molecule         The molecule.
     */
    FQMolecularHamiltonian(const Molecule& molecule);


    /*
     *  MARK: Access
     */

    /**
     *  @return The kinetic energy operator.
     */
    const KineticOperator& kinetic() const { return this->T; }

    /**
     *  @return The nuclear attraction operator.
     */
    const NuclearAttractionOperator& nuclearAttraction() const { return this->V; }

    /**
     *  @return The two-electron repulsion operator.
     */
    const CoulombRepulsionOperator& coulombRepulsion() const { return this->g; }
};


}  // namespace GQCP
