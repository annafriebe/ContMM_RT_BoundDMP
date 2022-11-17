/*!
 * @file    linear_qos_function.cpp
 * 
 * @brief   This class defines an implementation for the Quality of Service 
 *          functions.
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
#include "linear_qos_function.hpp"

namespace PrositCore {

  /// @brief Map linearly the quality.
  double LinearQoSFunction::evaluateQuality(double prob) {

    double qos = 0.0;

    /* The probability is smaller than the lower bound */
    if (prob <= pmin) {

      qos = offset;

    }
    else {

      /* The probability is greater than the upper bound */
      if (prob > pmax) {

        qos = (scale * (pmax - pmin)) + offset;

      }
      else {

        /* Return the quality based on a linear function */
        qos = (scale * (prob - pmin)) + offset;

      }

    }

    return qos;

  }

}
