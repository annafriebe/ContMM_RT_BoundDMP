/*!
 * @file    distribution.cpp
 * 
 * @brief   This class defines an implementation for probability distribution classes.
 * 
 * @author  Luigi Palopoli           <luigi.palopoli@unitn.it>
 *          Bernardo Villalba Frías  <b.r.villalba.frias@hva.nl>
 * 
 * @version 3.0
 * 
 * @date    30 November 2019
 * 
 * This program is a free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 */
#include "distribution.hpp"

using namespace std;

namespace PrositAux {

  /// @brief Load the file with the PMF of the distribution.
  int distribution::load(const string &filename) {

    /* Variables for the files */
    ifstream myfile;

    /* Variables for the data */
    int index;
    long double value;
  
    /* Open the data file */
    myfile.open(filename.c_str());
    if (myfile.fail()) {
      EXC_PRINT_2("ERROR: Unable to find the file: ", filename);
    }

    /* Check if the file was properly open */
    if (!myfile.is_open()) {
      EXC_PRINT_2("ERROR: Unable to open the file: ", filename);
    }

    /* Iterate over the file */
    while (myfile >> index >> value) {

      /* Set the element in the vector */
      set(index - min_val, value);
  
    }

    /* Close the file */
    myfile.close();

    return 0;

  }

  /// @brief Print out the distribution.
  void distribution::print() const {

    /* Iterate over the vector */
    for (uint64_t i = min_val; i <= max_val; i++) {

      /* Print out the vector */
      cout << i /* granularity*/ << " : " << setprecision(14) << elements(i - min_val) << endl;

    }

  }

  /// @brief Save the distribution.
  void distribution::dump(const std::string &filename) {

    /* Variables for the files */
    ofstream dumpFile;

    /* Open the data file */
    dumpFile.open(filename);

    /* Check if the file was properly open */
    if(dumpFile.is_open()) {

      /* Iterate over the vector */
      for (uint64_t i = min_val; i <= max_val; i++) {

        /* Save the vector */
        dumpFile << fixed << i << " " << setprecision(14) << elements(i - min_val) << endl;

      }

      /* Close the data file */
      dumpFile.close();

    }
    else {

      cerr << "WARNING: Unable to open the file: " << filename << endl;

    }

  }

}
