/*!
 * @file    qos_function.hpp
 * 
 * @brief   This class defines a header for the Quality of Service functions.
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
#ifndef QOS_FUNCTION_HPP
#define QOS_FUNCTION_HPP
#include <iomanip>
#include "../utils/exceptions.hpp"

namespace PrositCore {

  class QoSFunction {

    protected:

      double scale;     ///< Slope for the mapping.
      double pmin;      ///< Lower bound for the probability of deadline hit (quality offset below).
      double pmax;      ///< Upper bound for the probability of deadline hit (quality saturates above).

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for the QoS function.
      ///
      /// @param scaled is the slope for the linear function.
      /// @param pmind is the lower bound for the probability.
      /// @param pmaxd is the upper bound for the probability.
      /// @param offsetd is the minimum value for the quality.
      QoSFunction(double scaled,
          double pmind,
          double pmaxd) : 
          scale(scaled),
          pmin(pmind),
          pmax(pmaxd) {

        /* Check the correctness of the bounds */
        if ((pmin > pmax) || (pmin < 0) || (pmin > 1.0) || 
            (pmax < 0) || (pmax > 1.0)) {

          EXC_PRINT("ERROR: The probability limits are wrong.");

        }

        /* Check the slope of the function */
        if (scaled < 0) {

          EXC_PRINT("ERROR: The scaling constant is wrong.");

        }

      }

      /// @brief Destructor.
      ///
      /// This is the destructor for the QoS function.
      ~QoSFunction() { }

      /// @brief Map the quality.
      ///
      /// This function evaluates the quality of service based on the desired
      /// function.
      ///
      /// @param prob is the probability of respecting the deadline.
      ///
      /// @return the quality of the specified probability.
      virtual double evaluateQuality(double prob) = 0;

  };

}

#endif
