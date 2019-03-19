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
#ifndef GQCP_AOBASIS_HPP
#define GQCP_AOBASIS_HPP


#include "Basis/ShellSet.hpp"
#include "Operator/OneElectronOperator.hpp"



namespace GQCP {


/**
 *  A class that represents an atomic orbital basis, i.e. the collection of (scalar) atomic orbitals/basis functions
 */
class AOBasis {
private:
    ShellSet basisset;  // the underlying basisset that contains shells


public:
    // CONSTRUCTORS


    // PUBLIC METHODS - LIBINT INTEGRALS
    /**
     *  @return the matrix representation of the overlap operator in this AO basis
     */
    OneElectronOperator<double> calculateLibintOverlapIntegrals() const;


    // PUBLIC METHODS - LIBCINT INTEGRALS
};


}  // namespace GQCP


#endif  // GQCP_AOBASIS_HPP
