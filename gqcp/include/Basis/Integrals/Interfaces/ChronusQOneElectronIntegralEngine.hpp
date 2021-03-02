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


#include "Basis/Integrals/BaseOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/ChronusQ/engines.hpp"
#include "Basis/Integrals/Interfaces/ChronusQOneElectronIntegralBuffer.hpp"
#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/ScalarBasis/GIAOGTOShell.hpp"
#include "Utilities/aliases.hpp"


namespace GQCP {


/**
 *  A one-electron integral engine that uses ChronusQ as its back-end.
 * 
 *  @tparam _N              The number of components related to the operator.
 */
template <size_t _N>
class ChronusQOneElectronIntegralEngine:
    public BaseOneElectronIntegralEngine<GIAOGTOShell, _N, complex> {

public:
    // The scalar representation of a GIAO integral is a complex number.
    using IntegralScalar = complex;

    // The number of components related to the operator.
    static constexpr auto N = _N;

public:
    /*
     *  MARK: Constructors
     */

    /*
     *  MARK: Calculations
     */

    /**
     *  Calculate all the integrals over the given shells.
     * 
     *  @param shell1           The first shell.
     *  @param shell2           The second shell.
     * 
     *  @return A buffer containing the calculated integrals.
     * 
     *  @note This method is not marked const to allow the Engine's internals to be changed.
     */
    std::shared_ptr<BaseOneElectronIntegralBuffer<IntegralScalar, N>> calculate(const GIAOGTOShell& shell1, const GIAOGTOShell& shell2) override {

        // Prepare the ingredients for the ChronusQ `ComplexGIAOIntEngine`.
        auto libint_shell1 = LibintInterfacer::get().interface(shell1);
        auto libint_shell2 = LibintInterfacer::get().interface(shell2);
        libint2::ShellPair pair;
        pair.init(libint_shell1, libint_shell2, -1000);

        const auto& B = shell1.field().strength();
        std::array<double, 3> B_array {B(0), B(1), B(2)};

        const auto chronusq_buffer = ChronusQ::ComplexGIAOIntEngine::computeGIAOOverlapS(pair, libint_shell1, libint_shell2, B_array);

        return std::make_shared<ChronusQOneElectronIntegralBuffer<N>>(chronusq_buffer, shell1.numberOfBasisFunctions(), shell2.numberOfBasisFunctions());
    }
};


}  // namespace GQCP
