// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef SelectedFockSpace_hpp
#define SelectedFockSpace_hpp

#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/ONV.hpp"


namespace GQCP {


/**
 *  A class that represents a (spin-unresolved) Fock space that is flexible in the number of states (ONVs) that span it
 */
class SelectedFockspace : BaseFockSpace {
private:
    size_t M;  // number of spin orbitals
    size_t N;  // number of electrons

    std::vector<ONV> onvs;


public:

};



}  // namespace GQCP


#endif /* SelectedFockSpace_hpp */
