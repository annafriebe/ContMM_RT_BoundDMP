/*!
 * @file    companion_rr_probability_solver.hpp
 * 
 * @brief   This class defines the header for the Companion form solver for 
 *          resource reservation tasks described as a Quasi-Birth-Death 
 *          process.
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
#ifndef COMPANION_RR_PROBABILITY_SOLVER_HPP
#define COMPANION_RR_PROBABILITY_SOLVER_HPP

#include "../utils/auxiliary_functions.hpp"
#include "rr_probability_solver.hpp"

namespace PrositCore {

  class CompanionResourceReservationProbabilitySolver : public ResourceReservationProbabilitySolver {

    protected:

      Eigen::MatrixXd W;          ///< QBDP transition matrix in companion form.
      Eigen::MatrixXcd D;         ///< Power matrix W^j.
      uint64_t num_unst_eig;      ///< Number of unstable eigenvalues of W.
      Eigen::RowVectorXd Alfa;    ///< Vector with the alphas from the companion matrix.
      Eigen::RowVectorXd pi0;     ///< Stationary vector of the QBDP.
      //uint64_t max_deadline;      ///< Maximum deadline to compute.

      /// @brief Prepare the computation.
      ///
      /// This function generates the matrices used by the QBD methods.
      void preProcess();

      /// @brief Check the conditions to start.
      ///
      /// This function determines whether everything is ok to start the 
      /// pre-process phase.
      bool checkList();

      /// @brief Post processing after the solution.
      ///
      /// This function triggers the generation of the different probabilities.
      void postProcess() { }

      /// @brief Executes the iterations of the Companion form method.
      ///
      /// This function implements the Companion form method for discrete
      /// time QBDP.
      void applyAlgorithm();

      /// @brief Fills in the probability map.
      ///
      /// This function computes the probability map from the obtained 
      /// stationary vector pi0.
      void fillInProbabilityMap();

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for Companion form resource reservation 
      /// probability solvers.
      CompanionResourceReservationProbabilitySolver(/*uint64_t max_deadlined*/) : 
          W(),
          D(),
          num_unst_eig(0),
          Alfa(),
          pi0()/*,
          max_deadline(max_deadlined)*/ { }

      /// @brief Destructor.
      ///
      /// This is the destructor for Companion form resource reservation 
      /// probability solvers.
      virtual ~CompanionResourceReservationProbabilitySolver() { }

      std::complex<double> GammaFunctionsTotal(const std::complex<double> &b, uint64_t dim, const Eigen::RowVectorXd &Alpha);

  };

}

#endif
