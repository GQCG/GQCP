// This file is part of GQCG-template.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-template is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-template is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-template.  If not, see <http://www.gnu.org/licenses/>.


#include <boost/program_options.hpp>
#include <iomanip>
#include <iostream>
#include <string>

#include "temp.hpp"

namespace po = boost::program_options;


/**
 *  Main function for the executable
 */
int main (int argc, char** argv) {

    std::string filename_xyz;
    std::string basis_set;

    po::variables_map variables_map;

    /*
     *  Get the user's input
     */
    try {
        po::options_description desc ("Options");
        desc.add_options()
            ("help,h", "print help messages")
            ("filename_xyz,f", po::value<std::string>(&filename_xyz)->required(), "filename of the xyz file")
            ("basis_set,b", po::value<std::string>(&basis_set)->required(), "basis set name")
            ("charge,c", po::value<int>()->default_value(0), "charge")
            ("rdm,r", po::value<int>()->default_value(0), "reduced density matrices, 1 for 1-RDM, 2 for 2-RDM");

        po::store(po::parse_command_line(argc, argv, desc), variables_map);

        if (variables_map.count("help")) {
            std::cout << "DOCI" << std::endl << desc << std::endl;
            std::exit(0);
        }

        po::notify(variables_map);
    } catch (po::error& e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return 1;
    } catch(...) {
        std::cerr << "ERROR: you have not specified all arguments. Please use the -h flag for more information." << std::endl << std::endl;
    }

}  // main
