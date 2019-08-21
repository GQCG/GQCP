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
#include "Basis/Integrals/IntegralEngine.hpp"


namespace GQCP {


/*
 *  LIBINT - ONE-ELECTRON ENGINES
 */

/**
 *  @param op               the overlap operator
 *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
 *  @param max_l            the maximum angular momentum of Gaussian shell
 * 
 *  @return a one-electron integral engine that can calculate integrals over the overlap operator using the Libint integral library backend
 */
auto IntegralEngine::Libint(const OverlapOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<OverlapOperator::Components> {

    return LibintOneElectronIntegralEngine<OverlapOperator::Components>(op, max_nprim, max_l);
}


/**
 *  @param op               the kinetic operator
 *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
 *  @param max_l            the maximum angular momentum of Gaussian shell
 * 
 *  @return a one-electron integral engine that can calculate integrals over the kinetic operator using the Libint integral library backend
 */
auto IntegralEngine::Libint(const KineticOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<KineticOperator::Components> {

    return LibintOneElectronIntegralEngine<KineticOperator::Components>(op, max_nprim, max_l);
}


/**
 *  @param op               the nuclear attraction operator
 *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
 *  @param max_l            the maximum angular momentum of Gaussian shell
 * 
 *  @return a one-electron integral engine that can calculate integrals over the nuclear attraction operator using the Libint integral library backend
 */
auto IntegralEngine::Libint(const NuclearAttractionOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<NuclearAttractionOperator::Components> {

    return LibintOneElectronIntegralEngine<NuclearAttractionOperator::Components>(op, max_nprim, max_l);
}


/**
 *  @param op               the electronic electric dipole operator
 *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
 *  @param max_l            the maximum angular momentum of Gaussian shell
 * 
 *  @return a one-electron integral engine that can calculate integrals over the electronic electric dipole operator using the Libint integral library backend
 */
auto IntegralEngine::Libint(const ElectronicDipoleOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<ElectronicDipoleOperator::Components> {

    return LibintOneElectronIntegralEngine<ElectronicDipoleOperator::Components>(op, max_nprim, max_l);
}



/*
 *  LIBINT - TWO-ELECTRON ENGINES
 */

/**
 *  @param op               the Coulomb repulsion operator
 *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
 *  @param max_l            the maximum angular momentum of Gaussian shell
 * 
 *  @return a two-electron integral engine that can calculate integrals over the Coulomb repulsion operator using the Libint integral library backend
 */
auto IntegralEngine::Libint(const CoulombRepulsionOperator& op, const size_t max_nprim, const size_t max_l) -> LibintTwoElectronIntegralEngine<CoulombRepulsionOperator::Components> {

    return LibintTwoElectronIntegralEngine<CoulombRepulsionOperator::Components>(op, max_nprim, max_l);
}



/*
 *  LIBCINT - ONE-ELECTRON ENGINES
 */

/**
 *  @param op               the overlap operator
 *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
 * 
 *  @return a one-electron integral engine that can calculate integrals over the overlap operator using the Libcint integral library backend
 */
auto IntegralEngine::Libcint(const OverlapOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, OverlapOperator::Components, double> {

    return LibcintOneElectronIntegralEngine<GTOShell, OverlapOperator::Components, double>(op, shell_set);
}


/**
 *  @param op               the kinetic operator
 *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
 * 
 *  @return a one-electron integral engine that can calculate integrals over the kinetic operator using the Libcint integral library backend
 */
auto IntegralEngine::Libcint(const KineticOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, KineticOperator::Components, double> {

    return LibcintOneElectronIntegralEngine<GTOShell, KineticOperator::Components, double>(op, shell_set);
}


/**
 *  @param op               the nuclear attraction operator
 *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
 * 
 *  @return a one-electron integral engine that can calculate integrals over the nuclear attraction operator using the Libcint integral library backend
 */
auto IntegralEngine::Libcint(const NuclearAttractionOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, NuclearAttractionOperator::Components, double> {

    return LibcintOneElectronIntegralEngine<GTOShell, NuclearAttractionOperator::Components, double>(op, shell_set);
}


/**
 *  @param op               the electron electronic dipole operator
 * 
 *  @return a one-electron integral engine that can calculate integrals over the electronic dipole operator using the Libcint integral library backend
 */
auto IntegralEngine::Libcint(const ElectronicDipoleOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, ElectronicDipoleOperator::Components, double> {

    return LibcintOneElectronIntegralEngine<GTOShell, ElectronicDipoleOperator::Components, double>(op, shell_set);
}



/*
 *  LIBCINT - TWO-ELECTRON ENGINES
 */

/**
 *  @param op               the Coulomb repulsion operator
 *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
 * 
 *  @return a two-electron integral engine that can calculate integrals over the Coulomb repulsion operator using the Libcint integral library backend
 */
auto IntegralEngine::Libcint(const CoulombRepulsionOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintTwoElectronIntegralEngine<GTOShell, CoulombRepulsionOperator::Components, double> {

    return LibcintTwoElectronIntegralEngine<GTOShell, CoulombRepulsionOperator::Components, double>(op, shell_set);
}


}  // namespace GQCP
