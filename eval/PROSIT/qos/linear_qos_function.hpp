/*!
 * @file    linear_qos_function.hpp
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
#ifndef LINEAR_QOS_FUNCTION_HPP
#define LINEAR_QOS_FUNCTION_HPP

#include "qos_function.hpp"

namespace PrositCore {

  class LinearQoSFunction: public QoSFunction {

    private:

      double offset;    ///< Minimum value for the quality.

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for the QoS function.
      ///
      /// @param scaled is the slope for the linear function.
      /// @param pmind is the lower bound for the probability.
      /// @param pmaxd is the upper bound for the probability.
      /// @param offsetd is the minimum value for the quality.
      LinearQoSFunction(double scaled,
          double pmind,
          double pmaxd,
          double offsetd) : 
          QoSFunction(scaled, pmind, pmaxd),
          offset(offsetd) { }

      /// @brief Destructor.
      ///
      /// This is the destructor for the QoS function.
      ~LinearQoSFunction() { }

      /// @brief Map linearly the quality.
      ///
      /// This function evaluates the quality of service based on the desired
      /// linear function.
      ///
      /// @param prob is the probability of respecting the deadline.
      ///
      /// @return the quality of the specified probability.
      double evaluateQuality(double prob);

  };

}

#endif
