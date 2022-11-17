/*!
 * @file    analytic_rr_probability_solver.hpp
 * 
 * @brief   This class defines the header for the Analytic solver for resource
 *          reservation tasks described as a Quasi-Birth-Death process.
 * 
 *          [1] L. Palopoli, D. Fontanelli, L. Abeni and B. Villalba Frias. An
 *              analytical solution for probabilistic guarantees of reservation 
 *              based soft real-time systems. IEEE Transactions on Parallel and 
 *              Distributed Systems, vol. 27, no. 3, pp. 640–653, March 2016.
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
#ifndef ANALYTIC_RR_PROBABILITY_SOLVER_HPP
#define ANALYTIC_RR_PROBABILITY_SOLVER_HPP

#include "../utils/auxiliary_functions.hpp"
#include "rr_probability_solver.hpp"

namespace PrositCore {

  class AnalyticResourceReservationProbabilitySolver : public ResourceReservationProbabilitySolver {

    protected:

      Eigen::RowVectorXd Alfa;    ///< Vector with the alphas from the analytic form.
      Eigen::RowVectorXd pi0;     ///< Stationary vector of the QBDP.

      /// @brief Prepare the computation.
      ///
      /// This function generates the matrices required for the solution.
      void preProcess() { }

      /// @brief Post processing after the solution.
      ///
      /// This function checks that everything is ok for computation.
      void postProcess() { }

      /// @brief Fills in the probability map.
      ///
      /// This function computes the various probabilities after the "core" 
      /// problem has been solved.
      void fillInProbabilityMap();

      /// @brief Applies the algorithm.
      ///
      /// This function applies the selected algorithm to compute the
      /// probabilistic deadlines.
      void applyAlgorithm();

      void computeNextVector(Eigen::RowVectorXd &pi);

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for analytic resource reservation probability 
      /// solvers.
      ///
      /// @param Cd is the pointer to the distribution of the computation time.
      ///        The solver takes ownership of this pointer.
      /// @param Tsd is the reservation period.
      /// @param Qsd is the reservation budget.
      AnalyticResourceReservationProbabilitySolver() : 
          Alfa(),
          pi0() { }

      /// @brief Destructor.
      ///
      /// This is the destructor for analytic resource reservation probability 
      /// solvers.
      virtual ~AnalyticResourceReservationProbabilitySolver() { }

  };

}

#endif
