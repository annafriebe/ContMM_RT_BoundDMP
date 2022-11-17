/*!
 * @file    cr_qbd_rr_probability_solver.hpp
 * 
 * @brief   This class defines the header for the Cyclic Reduction [1] solver
 *          for resource reservation tasks based on algorithms for numeric 
 *          solution of Quasi-Birth-Death processes.
 * 
 *          [1] B. Meing and D. Bini. Improved cyclic reduction for solving 
 *              queueing problems. Numerical Algorithms, 15:57--74, 1997.
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
#ifndef CR_QBD_RR_PROBABILITY_SOLVER_HPP
#define CR_QBD_RR_PROBABILITY_SOLVER_HPP

#include "qbd_rr_probability_solver.hpp"

namespace PrositCore {

  class CRQBDResourceReservationProbabilitySolver : public QBDResourceReservationProbabilitySolver {

    protected:

      bool shift_flag;       ///< Shift technique to accelerate convergence.
      int max_iter;          ///< Maximum number of iteration of the algorithm.
      Eigen::MatrixXd G;     ///< Minimal nonnegative solution to G = A0 + A1 G + A2 G^2.
      Eigen::MatrixXd U;     ///< Minimal nonnegative solution to U = A1 + A2 (I-U)^(-1) A0.

      /// @brief Executes the iterations of the Cyclic Reduction method.
      ///
      /// This function implements the Cyclic Reduction method for discrete
      /// time QBDP. The iteration produces a matrix R that has to be 
      /// post-processed.
      void applyAlgorithm();

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for Cyclic Reduction QBD resource reservation
      /// probability solvers.
      ///
      /// @param shift_d is the shift parameter.
      /// @param max_iter_d is the number of iterations.
      CRQBDResourceReservationProbabilitySolver(bool shift_d, int max_iter_d) : 
          QBDResourceReservationProbabilitySolver(), 
          shift_flag(shift_d), 
          max_iter(max_iter_d), 
          G(), 
          U() { }

      /// @brief Destructor.
      ///
      /// This is the destructor for Cyclic Reduction QBD resource reservation 
      /// probability solvers.
      virtual ~CRQBDResourceReservationProbabilitySolver() { }

      /// @brief Sets the drift parameter.
      ///
      /// This function sets to true the shift flag.
      void setShift() { 
        shift_flag = true;
      }

      /// @brief Sets the maximum number of iterations.
      ///
      /// This function sets the maximum number of iterations of the algorithm.
      ///
      /// @param max_iter_d is the maximum number of iterations.
      void setMaxIter(int max_iter_d) { 
        max_iter = max_iter_d; 
      }

  };

}

#endif
