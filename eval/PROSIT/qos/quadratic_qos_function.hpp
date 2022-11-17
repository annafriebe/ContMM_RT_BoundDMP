/*!
 * @file    quadratic_qos_function.hpp
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
#ifndef QUADRATIC_QOS_FUNCTION_HPP
#define QUADRATIC_QOS_FUNCTION_HPP

#include "qos_function.hpp"

namespace PrositCore {

  class QuadraticQoSFunction: public QoSFunction {

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for the QoS function.
      ///
      /// @param scaled is the slope for the linear function.
      /// @param pmind is the lower bound for the probability.
      /// @param pmaxd is the upper bound for the probability.
      QuadraticQoSFunction(double scaled,
          double pmind,
          double pmaxd) : 
          QoSFunction(scaled, pmind, pmaxd) { }

      /// @brief Destructor.
      ///
      /// This is the destructor for the QoS function.
      ~QuadraticQoSFunction() { }

      /// @brief Map quadratically the quality.
      ///
      /// This function evaluates the quality of service based on the desired
      /// quadratic function.
      ///
      /// @param prob is the probability of respecting the deadline.
      ///
      /// @return the quality of the specified probability.
      double evaluateQuality(double prob);

  };

}

#endif
