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


#include "Operator/FirstQuantized/AngularMomentumOperator.hpp"
#include "Operator/FirstQuantized/BaseFQOperator.hpp"
#include "Operator/FirstQuantized/BaseReferenceDependentOperator.hpp"
#include "Physical/HomogeneousMagneticField.hpp"


namespace GQCP {


/**
 *  The paramagnetic operator, i.e. the part of the scalar kinetic energy operator that is linear in the magnetic field.
 */
class ParamagneticOperator:
    public BaseScalarFQOneElectronOperator<complex>,
    public BaseReferenceDependentOperator {
private:
    // The external, homogeneous magnetic field.
    HomogeneousMagneticField B;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a `ParamagneticOperator` from its underlying homogeneous magnetic field and a reference point.
     * 
     *  @param B                    The external, homogeneous magnetic field.
     *  @param reference            The reference point about which the paramagnetic operator is calculated.
     */
    ParamagneticOperator(const HomogeneousMagneticField& B, const Vector<double, 3>& reference);

    /**
     *  Construct a `ParamagneticOperator` from its underlying homogeneous magnetic field.
     * 
     *  @param B                    The external, homogeneous magnetic field.
     * 
     *  @note The reference point is chosen to be the gauge origin of the magnetic field.
     */
    ParamagneticOperator(const HomogeneousMagneticField& B);


    /*
     *  MARK: Access
     */

    /**
     *  @return The external, homogeneous magnetic field.
     */
    const HomogeneousMagneticField& magneticField() const { return this->B; }


    /*
     *  MARK: First-quantized operators
     */

    /**
     *  @return The angular momentum operator that is related to this paramagnetic operator.
     */
    AngularMomentumOperator angularMomentum() const { return AngularMomentumOperator(this->reference()); }
};


}  // namespace GQCP
