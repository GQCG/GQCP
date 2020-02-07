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
#pragma once


#include "ONVBasis/SelectedONVBasis.hpp"
#include "QCModel/CI/LinearExpansion.hpp"


namespace GQCP {


/**
 *  A class that reads and stores a 'selected' wave function expansion
 */
class LinearExpansionReader {
private:
    SelectedONVBasis fock_space;
    VectorX<double> coefficients;
    LinearExpansion wave_function;


public:
    /**
     *  @param GAMESS_filename      the name of the GAMESS file that contains the 'selected' wave function expansion
     */
    explicit LinearExpansionReader(const std::string& GAMESS_filename);


    // GETTERS
    const SelectedONVBasis& get_fock_space() const { return this->fock_space; }
    const VectorX<double>& get_coefficients() const { return this->coefficients; }
    const LinearExpansion& get_wave_function() const { return this->wave_function; }
};



}  // namespace GQCP
