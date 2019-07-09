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
#ifndef GQCP_LIBINTONEELECTRONINTEGRALENGINE_HPP
#define GQCP_LIBINTONEELECTRONINTEGRALENGINE_HPP


#include "Basis/BaseOneElectronIntegralEngine.hpp"

#include "Basis/GTOShell.hpp"
#include "Basis/LibintInterfacer.hpp"
#include "Basis/LibintOneElectronIntegralBuffer.hpp"
#include "Operator/FirstQuantized/Operator.hpp"



namespace GQCP {


/**
 *  A one-electron integral engine that uses libint as its backend
 * 
 *  @tparam _N                  the number of components the operator has
 */
template <size_t _N>
class LibintOneElectronIntegralEngine : public BaseOneElectronIntegralEngine<GTOShell, _N, double> {
public:
    using IntegralScalar = double;  // the scalar representation of an integral for libint is always a real number
    static constexpr auto N = _N;  // the number of components the operator has


private:
    libint2::Engine libint2_engine;


    // Parameters to give to the buffer
    size_t component_offset = 0;  // the number of libint components that should be skipped during access of calculated values (libint2::Operator::emultipole1 has 4 libint2 components, but in reality there should only be 1)
    double scaling_factor = 1.0;  // a factor that is multiplied to all of the calculated integrals


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param op               the overlap operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     */
    LibintOneElectronIntegralEngine(const OverlapOperator& op, const size_t max_nprim, const size_t max_l) :
        libint2_engine (LibintInterfacer::get().createEngine(op, max_nprim, max_l))
    {}

    /**
     *  @param op               the kinetic operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     */
    LibintOneElectronIntegralEngine(const KineticOperator& op, const size_t max_nprim, const size_t max_l) :
        libint2_engine (LibintInterfacer::get().createEngine(op, max_nprim, max_l))
    {}

    /**
     *  @param op               the nuclear attraction operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     */
    LibintOneElectronIntegralEngine(const NuclearAttractionOperator& op, const size_t max_nprim, const size_t max_l) :
        libint2_engine (LibintInterfacer::get().createEngine(op, max_nprim, max_l))
    {
        auto libint_atoms = LibintInterfacer::get().interface(op.nuclearFramework().nucleiAsVector());
        this->libint2_engine.set_params(libint2::make_point_charges(libint_atoms));
    }

    /**
     *  @param op               the electronic electric dipole operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     */
    LibintOneElectronIntegralEngine(const ElectronicDipoleOperator& op, const size_t max_nprim, const size_t max_l) :
        libint2_engine (LibintInterfacer::get().createEngine(op, max_nprim, max_l)),
        component_offset (1),  // emultipole1 has [overlap, x, y, z], we don't need the overlap
        scaling_factor (-1.0)  // apply the minus sign which comes from the charge of the electrons -e
    {
        std::array<double, 3> libint2_origin_array {op.origin().x(), op.origin().y(), op.origin().z()};
        this->libint2_engine.set_params(libint2_origin_array);
    }



    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @param shell1           the first shell
     *  @param shell2           the second shell
     */
    std::shared_ptr<BaseOneElectronIntegralBuffer<IntegralScalar, N>> calculate(const GTOShell& shell1, const GTOShell& shell2) override {

        const auto libint_shell1 = LibintInterfacer::get().interface(shell1);
        const auto libint_shell2 = LibintInterfacer::get().interface(shell2);

        const auto& libint2_buffer = this->libint2_engine.results();
        this->libint2_engine.compute(libint_shell1, libint_shell2);
        return std::make_shared<LibintOneElectronIntegralBuffer<N>>(libint2_buffer, shell1.numberOfBasisFunctions(), shell2.numberOfBasisFunctions(), this->component_offset, this->scaling_factor);
    }
};


}  // namespace GQCP



#endif  // GQCP_LIBINTONEELECTRONINTEGRALENGINE_HPP
