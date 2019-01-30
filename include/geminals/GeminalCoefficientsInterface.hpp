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
#ifndef GeminalCoefficientsInterface_hpp
#define GeminalCoefficientsInterface_hpp


#include "WaveFunction/WaveFunction.hpp"
#include "FockSpace/FockSpace.hpp"


namespace GQCP {


class GeminalCoefficientsInterface {
public:
    // DESTRUCTOR
    virtual ~GeminalCoefficientsInterface();


    // PUBLIC METHODS
    /**
     *  @param onv      the ONV that is being projected on
     *
     *  @return the overlap of the APIG-like wave function with the given on, i.e. the projection of the APIG wave function onto that ONV
     */
    virtual double overlap(const ONV& onv) const = 0;

    /**
     *  @param fock_space       the seniority-zero Fock space the wave function should live in
     *
     *  @return the wave function expansion corresponding to the geminal coefficients
     */
    WaveFunction toWaveFunction(const FockSpace& fock_space) const;

};

}  // namespace GQCP



#endif /* GeminalCoefficientsInterface_hpp */
