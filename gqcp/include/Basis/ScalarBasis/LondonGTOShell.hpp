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


#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Mathematical/Functions/LondonCartesianGTO.hpp"
#include "Physical/HomogeneousMagneticField.hpp"


namespace GQCP {


/**
 *  A GTO shell with London-type modifications.
 */
class LondonGTOShell {
public:
    // The type of primitive that underlies this shell.
    using Primitive = LondonCartesianGTO;

    // The type of basis function that this shell can produce.
    using BasisFunction = LinearCombination<double, Primitive>;


private:
    // The specificiation of the magnetic field for this GIAO.
    HomogeneousMagneticField m_field;

    // The field-free GTO shell.
    GTOShell gto_shell;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a gauge-including shell from a field-free shell and a magnetic field.
     * 
     *  @param gto_shell        The field-free GTO shell.
     *  @param field            The specificiation of the magnetic field for this GIAO.
     */
    LondonGTOShell(const GTOShell& shell, const HomogeneousMagneticField& field);


    /*
     *  MARK: Access
     */

    /**
     *  @return The specificiation of the magnetic field for this GIAO.
     */
    const HomogeneousMagneticField& magneticField() const { return this->m_field; }


    /**
     *  @return The field-free GTO shell.
     */
    const GTOShell& gtoShell() const { return this->gto_shell; }


    /*
     *  MARK: Normalization
     */

    /**
     *  Embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing.
     *
     *  @note The normalization factor that is embedded, corresponds to the spherical (or axis-aligned Cartesian) GTO.
     */
    void embedNormalizationFactorsOfPrimitives() { this->gto_shell.embedNormalizationFactorsOfPrimitives(); }


    /*
     *  MARK: Basis functions
     */

    /**
     *  @return The number of basis functions that this shell contains.
     */
    size_t numberOfBasisFunctions() const { return this->gtoShell().numberOfBasisFunctions(); }

    /**
     *  Construct all basis functions contained in this shell.
     * 
     *  @return The basis functions that correspond to this shell.
     * 
     *  @note The basis functions are ordered lexicographically. This means x < y < z.
     */
    std::vector<BasisFunction> basisFunctions() const;
};


}  // namespace GQCP