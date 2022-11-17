/*!
 * @file    lt_qbd_rr_probability_solver.hpp
 * 
 * @brief   This class defines the header for the Latouche Reduction [1] solver
 *          for resource reservation tasks based on algorithms for numeric 
 *          solution of Quasi-Birth-Death processes.
 * 
 *          [1] G. Latouche and V. Ramaswami. A logarithmic reduction algorithm 
 *              for Quasi-Birth-and-Death processes. J. Appl. Probab., 30:650-674, 
 *              1993.
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
#ifndef LT_QBD_RR_PROBABILITY_SOLVER_HPP
#define LT_QBD_RR_PROBABILITY_SOLVER_HPP

#include "qbd_rr_probability_solver.hpp"

namespace PrositCore {

  class LatoucheQBDResourceReservationProbabilitySolver : public QBDResourceReservationProbabilitySolver {

    protected:

      double epsilon;    ///< Threshold for the result of two iterations.
      int max_iter;      ///< Maximum number of iteration of the algorithm.

      /// @brief Executes the iterations of the Latouche method.
      ///
      /// This function implements the Cyclic Reduction method for discrete
      /// time QBDP. The iteration produces a matrix R that has to be 
      /// post-processed.
      void applyAlgorithm();

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for Latouche Reduction QBD resource reservation
      /// probability solvers.
      ///
      /// @param epsilon_d is the desired value for the epsilon parameter.
      /// @param max_iter_d is the number of iterations.
      LatoucheQBDResourceReservationProbabilitySolver(double epsilon_d, unsigned int max_iter_d) : 
          QBDResourceReservationProbabilitySolver(), 
          epsilon(epsilon_d), 
          max_iter(max_iter_d) {

        /* Check the correctness of the epsilon */
        if (epsilon < 0) {
          EXC_PRINT("ERROR: The epsilon parameter has to be non negative.");
        }

      }

      /// @brief Destructor.
      ///
      /// This is the destructor for Latouche Reduction QBD resource reservation
      /// probability solvers.
      virtual ~LatoucheQBDResourceReservationProbabilitySolver() { }

      /// @brief Set the epsilon to a desired value.
      ///
      /// This function sets the epsilon value for the result of two consecutive
      /// iterations.
      ///
      /// @param epsilon_d is the desired value of epsilon.
      void setEpsilon(double epsilon_d) {

        /* Check the correctness of the epsilon */
        if (epsilon_d < 0) {
          EXC_PRINT("ERROR: The epsilon parameter has to be non negative.");
        }

        epsilon = epsilon_d;

      }

      /// @brief Sets the maximum number of iterations.
      ///
      /// This function sets the maximum number of iterations of the algorithm.
      ///
      /// @param max_iter_d is the maximum number of iterations.
      void setMaxIter(unsigned int max_iter_d) {
        max_iter = max_iter_d;
      }

  };

}

#endif
