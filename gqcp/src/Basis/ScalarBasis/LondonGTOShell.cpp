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

#include "Basis/ScalarBasis/LondonGTOShell.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  Construct a gauge-including shell from a field-free shell and a magnetic field.
 * 
 *  @param gto_shell        The field-free GTO shell.
 *  @param field            The specificiation of the magnetic field for this GIAO.
 */
LondonGTOShell::LondonGTOShell(const GTOShell& gto_shell, const HomogeneousMagneticField& field) :
    gto_shell {gto_shell},
    m_field {field} {}


/*
 *  MARK: Basis functions
 */

/**
 *  Construct all basis functions contained in this shell.
 * 
 *  @return The basis functions that correspond to this shell.
 * 
 *  @note The basis functions are ordered lexicographically. This means x < y < z.
 */
std::vector<LondonGTOShell::BasisFunction> LondonGTOShell::basisFunctions() const {

    // Every primitive in the basis functions of the Cartesian GTO shell needs to acquire the London modification.
    const auto gto_basis_functions = this->gtoShell().basisFunctions();

    std::vector<LondonGTOShell::BasisFunction> london_basis_functions {};
    london_basis_functions.reserve(this->numberOfBasisFunctions());
    for (const auto& bf : gto_basis_functions) {
        const auto& coefficients = bf.coefficients();
        const auto& gto_primitives = bf.functions();

        LondonGTOShell::BasisFunction london_basis_function {};
        for (size_t c = 0; c < bf.length(); c++) {
            const auto& coefficient = coefficients[c];
            const auto& gto_primitive = gto_primitives[c];

            const LondonCartesianGTO london_primitive {this->magneticField(), gto_primitive};
            london_basis_function.append(coefficient, london_primitive);
        }

        london_basis_functions.push_back(london_basis_function);
    }

    return london_basis_functions;
}


}  // namespace GQCP
