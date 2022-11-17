/*!
 * @file    pmf.cpp
 * 
 * @brief   This class defines an implementation for probability mass function 
 *          classes.
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
#include "pmf.hpp"

using namespace std;

namespace PrositAux {

  /// @brief Obtain the sum of the PMF.
  double pmf::sum() const {

    /* Variable to accumulate the sum */
    double p;

    /* Initially the sum contains the tail of the distribution */
    p = tail;

    /* Add the sum of the vector */
    p += elements.sum();

    return p;

  }

  /// @brief Obtain the mean of the PMF.
  double pmf::avg() const {

    /* Variable for the mean */
    double avg = 0;

    /* Iterate over the vector */
    for (uint64_t i = min_val; i <= max_val; i++) {

      /* Accumulate the mean */
      avg += (i * elements(i - min_val));

    }

    return avg;

  }

  /// @brief Obtain the standard deviation of the PMF.
  double pmf::std() const {

    /* Obtain the variance */
    double var = pmf::var();

    /* Return the standard deviation */
    return sqrt(var);

  }

  /// @brief Obtain the variance of the PMF.
  double pmf::var() const {

    /* Variables for the variance */
    double var = 0;
    double avg = pmf::avg();

    /* Iterate over the vector */
    for (uint64_t i = min_val; i <= max_val; i++) {

      /* Accumulate the variance */
      var += (pow((i - avg), 2) * elements(i - min_val));

    }

    return var;

  }

  /// @brief Set the value for an index in the distribution.
  int pmf::set(uint64_t idx, double val) {

    /* Check that the index is in range */
    if (idx > size) {
      EXC_PRINT("ERROR: Access out of range while setting.");
    }

    /* Set the element */
    elements(idx) = val;

    return 0;

  }

  /// @brief Get the value of the distribution.
  double pmf::get(uint64_t idx) const {

    /* Check that the index is in range */
    if (idx > (size - 1)) {
      EXC_PRINT("ERROR: Access out of range while getting.");
    }

    /* Return the value */
    return elements(idx);

  }

  /// @brief Check the PMF of the distribution.
  pmf::ERR_CODES pmf::check() const {

    /* Obtain the sum of the PMF */
    double p = sum();

    /* The PMF sums less than one */
    if (p < (1.0 - epsilon)) {
      cerr << "WARNING: The sum of the PMF is: " << p << endl;
      return PMF_SMALLER_ONE;
    }

    /* The PMF sums more than one */
    if (p > 1.0 + epsilon) {
      cerr << "WARNING: The sum of the PMF is: " << p << endl;
      return PMF_GREATER_ONE;
    }

    /* The PMF is correct */
    return PMF_OK;

  }

  /// @brief Load the file with the PMF of the distribution.
  int pmf::load(const string &filename) {

    /* Load the PMF */
    int j = distribution::load(filename);

    /* Check possible errors with the PMF */
    ERR_CODES e = check();

    /* Present the corresponding error messages */
    switch (e) {

      /* Case of PMF smaller than one */
      case PMF_SMALLER_ONE:
        EXC_PRINT_2("ERROR: The sum of the PMF is smaller than one, Check the file: ", filename);
        break;

      /* Case of PMF greater than one */
      case PMF_GREATER_ONE:
        EXC_PRINT_2("ERROR: The sum of the PMF is greater than one, Check the file: ", filename);
        break;

      default:
        ;

    }

    return j;

  }

  /// @brief Resample the PMF to a given granularity
  unique_ptr<pmf> pmf::resample(uint64_t delta, bool growing) {

    /* Initialize the variables */
    uint64_t lim = 0;
    double acum_prob = 0.0;

    /* Quantise to the upper limit */
    if (growing) {

      /* Create the proper new PMF */
      unique_ptr<pmf> clone(new pmf(ceil((double)this->getMin() / (double)delta),
          ceil((double)this->getMax() / (double)delta)/*, delta*/, this->epsilon));

      /* Obtain the first bin */
      lim = ceil((double)min_val / (double)delta) * delta;

      /* Iterate over the rows of the PMF */
      for (uint64_t i = min_val; i <= max_val; i++) {

        if (i <= lim) {
          acum_prob += get(i - min_val);
        }
        else {
          clone->set((lim / delta) - clone->getMin(), acum_prob);
          lim += delta;

          while (i > lim) {
            clone->set((lim / delta) - clone->getMin(), 0);
            lim += delta;
          }

          acum_prob = get(i - min_val);

        }

      }

      if (acum_prob != 0)
        clone->set((lim / delta) - clone->getMin(), acum_prob);

      return clone;

    }
    else {

      /* Create the proper new PMF */
      unique_ptr<pmf> clone(new pmf(floor((double)this->getMin() / (double)delta), 
          floor((double)this->getMax() / (double)delta)/*, delta*/, this->epsilon));

      /* Obtain the first bin */
      lim = floor((double)min_val / (double)delta) * delta;

      /* Iterate over the rows of the pmf file */
      for (uint64_t i = min_val; i <= max_val; i++) {

        if (i <= lim + delta - 1) {
          acum_prob += get(i - min_val);
        }
        else {
          clone->set((lim / delta) - clone->getMin(), acum_prob);
          lim += delta;

          while (i > lim + delta - 1) {
            clone->set((lim / delta) - clone->getMin(), 0);
            lim += delta;
          }

          acum_prob = get(i - min_val);

        }

      }

      if (acum_prob != 0)
        clone->set((lim / delta) - clone->getMin(), acum_prob);

      return clone;

    }

  }

}
