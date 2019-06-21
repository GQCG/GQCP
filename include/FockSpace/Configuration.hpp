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
#ifndef GQCP_CONFIGURATION_HPP
#define GQCP_CONFIGURATION_HPP


#include "ONV.hpp"


namespace GQCP {


<<<<<<< HEAD:include/FockSpace/Configuration.hpp
/**
 *  A struct that represents the combination of an alpha and beta ONV
 */
struct Configuration {
    ONV onv_alpha;
    ONV onv_beta;
=======
class RHFElectricalResponseSolver : BaseElectricalResponseSolver {
private:
    size_t N_P;  // the number of electron pairs


public:
    // CONSTRUCTORS

    /**
     *  @param N_P          the number of electron pairs
     */
    RHFElectricalResponseSolver(const size_t N_P);


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param ham_par                  the Hamiltonian parameters
     * 
     *  @return the parameter response constant (k_p), i.e. the second-order parameter partial derivative of the RHF energy function
     */
    SquareMatrix<double> calculateParameterResponseConstant(const HamiltonianParameters<double>& ham_par) override;

    /**
     *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
     * 
     *  @return the parameter response force (F_p), i.e. the first-order parameter partial derivative of the perturbation derivative of the RHF energy function
     */
    Matrix<double, Dynamic, 3> calculateParameterResponseForce(const std::array<OneElectronOperator<double>, 3>& dipole_integrals) override;
>>>>>>> Implement RHF linear response equations:include/Properties/RHFElectricalResponseSolver.hpp
};


}  // namespace GQCP


#endif  // GQCP_CONFIGURATION_HPP
