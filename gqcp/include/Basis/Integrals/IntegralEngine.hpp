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


#include "Basis/Integrals/Interfaces/LibcintOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibcintTwoElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibintOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibintTwoElectronIntegralEngine.hpp"
#include "Basis/Integrals/OneElectronIntegralEngine.hpp"
#include "Basis/Integrals/PrimitiveAngularMomentumIntegralEngine.hpp"
#include "Basis/Integrals/PrimitiveDipoleIntegralEngine.hpp"
#include "Basis/Integrals/PrimitiveKineticEnergyIntegralEngine.hpp"
#include "Basis/Integrals/PrimitiveLinearMomentumIntegralEngine.hpp"
#include "Basis/Integrals/PrimitiveOverlapIntegralEngine.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Utilities/aliases.hpp"


namespace GQCP {


/**
 *  A class that creates integral engines, like a factory class.
 */
class IntegralEngine {
public:
    /*
     *  GQCP ("In-house")
     */

    /**
     *  @param op               the angular momentum operator
     * 
     *  @return a one-electron integral engine that can calculate integrals over the angular momentum operator
     */
    static OneElectronIntegralEngine<PrimitiveAngularMomentumIntegralEngine> InHouse(const AngularMomentumOperator& op);

    /**
     *  @param op               the electronic dipole operator
     * 
     *  @return a one-electron integral engine that can calculate integrals over the electronic dipole operator
     */
    static OneElectronIntegralEngine<PrimitiveDipoleIntegralEngine> InHouse(const ElectronicDipoleOperator& op);

    /**
     *  @param op               the kinetic energy operator
     * 
     *  @return a one-electron integral engine that can calculate integrals over the kinetic energy operator
     */
    static OneElectronIntegralEngine<PrimitiveKineticEnergyIntegralEngine> InHouse(const KineticOperator& op);

    /**
     *  @param op               the linear momentum operator
     * 
     *  @return a one-electron integral engine that can calculate integrals over the linear momentum operator
     */
    static OneElectronIntegralEngine<PrimitiveLinearMomentumIntegralEngine> InHouse(const LinearMomentumOperator& op);

    /**
     *  @param op               the overlap operator
     * 
     *  @return a one-electron integral engine that can calculate integrals over the overlap operator
     */
    static OneElectronIntegralEngine<PrimitiveOverlapIntegralEngine> InHouse(const OverlapOperator& op);


    /*
     *  LIBINT
     */

    /**
     *  @param op               the Coulomb repulsion operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return a two-electron integral engine that can calculate integrals over the Coulomb repulsion operator using the Libint integral library backend
     */
    static auto Libint(const CoulombRepulsionOperator& op, const size_t max_nprim, const size_t max_l) -> LibintTwoElectronIntegralEngine<CoulombRepulsionOperator::Components>;

    /**
     *  @param op               the electronic electric dipole operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return a one-electron integral engine that can calculate integrals over the electronic electric dipole operator using the Libint integral library backend
     */
    static auto Libint(const ElectronicDipoleOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<ElectronicDipoleOperator::Components>;

    /**
     *  @param op               the nuclear attraction operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return a one-electron integral engine that can calculate integrals over the nuclear attraction operator using the Libint integral library backend
     */
    static auto Libint(const NuclearAttractionOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<NuclearAttractionOperator::Components>;

    /**
     *  @param op               the kinetic operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return a one-electron integral engine that can calculate integrals over the kinetic operator using the Libint integral library backend
     */
    static auto Libint(const KineticOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<KineticOperator::Components>;

    /**
     *  @param op               the overlap operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return a one-electron integral engine that can calculate integrals over the overlap operator using the Libint integral library backend
     */
    static auto Libint(const OverlapOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<OverlapOperator::Components>;


    /*
     *  LIBCINT
     */

    /**
     *  @param op               the Coulomb repulsion operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     * 
     *  @return a two-electron integral engine that can calculate integrals over the Coulomb repulsion operator using the Libcint integral library backend
     */
    static auto Libcint(const CoulombRepulsionOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintTwoElectronIntegralEngine<GTOShell, CoulombRepulsionOperator::Components, double>;

    /**
     *  @param op               the electron electronic dipole operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     * 
     *  @return a one-electron integral engine that can calculate integrals over the electronic dipole operator using the Libcint integral library backend
     */
    static auto Libcint(const ElectronicDipoleOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, ElectronicDipoleOperator::Components, double>;

    /**
     *  @param op               the nuclear attraction operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     * 
     *  @return a one-electron integral engine that can calculate integrals over the nuclear attraction operator using the Libcint integral library backend
     */
    static auto Libcint(const NuclearAttractionOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, NuclearAttractionOperator::Components, double>;

    /**
     *  @param op               the kinetic operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     * 
     *  @return a one-electron integral engine that can calculate integrals over the kinetic operator using the Libcint integral library backend
     */
    static auto Libcint(const KineticOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, KineticOperator::Components, double>;

    /**
     *  @param op               the overlap operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     * 
     *  @return a one-electron integral engine that can calculate integrals over the overlap operator using the Libcint integral library backend
     */
    static auto Libcint(const OverlapOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, OverlapOperator::Components, double>;
};


}  // namespace GQCP
