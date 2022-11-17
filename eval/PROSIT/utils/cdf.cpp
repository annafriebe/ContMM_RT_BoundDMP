/*!
 * @file    cdf.cpp
 * 
 * @brief   This class defines an implementation for cumulative distribution 
 *          function classes.
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
#include "cdf.hpp"

using namespace std;

namespace PrositAux {

  /// @brief Constructor.
  cdf::cdf(int cmin, int cmax/*, int delta*/, double eps) :
      distribution(cmin, cmax/*, delta*/, eps), just_created(true) {

    /* Initialize the vector with the distribution */
    elements.setZero();

  }

  /// @brief Get the value of the distribution.
  double cdf::get(uint64_t idx) const {

    /* Checks that the cdf has data */
    if (just_created) {
      EXC_PRINT("ERROR: Function called on an empty CDF.");
    }

    /* The given index is greater than the maximum */
    if (idx > size) {

      /* The CDF is correct */
      if (elements(max_val - min_val) > 1.0 - epsilon)
        return 1.0;
      else
        EXC_PRINT("ERROR: The CDF is ill formed.");

    }

    /* Return the value */
    return elements(idx);

  }

  /// @brief Load the file with the CDF of the distribution.
  int cdf::load(const string &filename) {

    /* Load the CDF */
    int j = distribution::load(filename);

    /* Check possible errors with the PMF */
    if (check() != CDF_OK)
      EXC_PRINT_2("ERROR: Ill formed cdf", filename);

    return j;

  }

  /// @brief Check the CDF of the distribution.
  cdf::ERR_CODES cdf::check() const {

    /* Checks that the cdf has data */
    if (just_created)
      EXC_PRINT("ERROR: Function called on an empty CDF.");

    /* The CDF does not reach 1 */
    if ((elements(max_val - min_val) > 1 + epsilon) || (elements(max_val - min_val) < 1 - epsilon)) {
      cerr << "WARNING: The maximum is incorrect" << endl;
      return CDF_BAD_MAX;
    }

    /* The CDF is not increasing */
    for (uint64_t i = 1; i < size; i++) {
      if (get(i - 1) > get(i))
        return CDF_NON_INCREASING;
    }

    /* The CDF is correct */
    return CDF_OK;

  }

  /// @brief Set the value for an index in the distribution.
  int cdf::set(uint64_t idx, double val) {

    /* Check that the index is in range */
    if (idx > size)
      EXC_PRINT("ERROR: Insert out of range.");

    /* Set the proper flag */
    if (just_created)
      just_created = false;

    /* Set the element */
    elements(idx) = val;

    return 0;

  }

}
