/*!
 * @file    probability_solver.hpp
 * 
 * @brief   This class contains the abstract definitions for probability 
 *          solvers. In essence, a probability solver wraps the different 
 *          algorithms used to compute probabilistic deadlines.
 *          Specific interfaces are also provided for resource reservations 
 *          and fixed priority.
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
#ifndef PROBABILITY_SOLVER_HPP
#define PROBABILITY_SOLVER_HPP

#include "../tasks/generic_task_descriptor.hpp"

namespace PrositCore {

  class ProbabilitySolver {

    protected:

      bool pre_process_done;    ///< The pre-process has been performed.
      bool solved;              ///< The solver has already been called.
      bool post_process_done;   ///< The post-process has been performed.
      bool verbose_flag;        ///< Flag to print out online information.
      bool linked_flag;         ///< The solver has been linked with one task.
  
      /// @brief Prepare the computation.
      ///
      /// This function prepares the preliminary elements to be used by the 
      /// solver.
      virtual void preProcess() = 0;

      /// @brief Applies the algorithm.
      ///
      /// This function applies the selected algorithm to compute the
      /// probabilistic deadlines.
      virtual void applyAlgorithm() = 0;

      /// @brief Post processing after the solution.
      ///
      /// This function completes the computation of the algorithm.
      virtual void postProcess() = 0;

      /// @brief Fills in the probability map.
      ///
      /// This function fills the probability map with the obtained 
      /// probabilities.
      virtual void fillInProbabilityMap() = 0;

      /// @brief Check the conditions to start.
      ///
      /// This function determines whether everything is ok to start the 
      /// pre-process phase.
      virtual bool checkList() = 0;

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for probability solvers.
      ProbabilitySolver() : 
          pre_process_done(false), 
          solved(false), 
          post_process_done(false), 
          verbose_flag(false), 
          linked_flag(false) { }

      /// @brief Destructor.
      ///
      /// This is the destructor of the probability solvers. It is defined
      /// virtual to enable subclassing.
      virtual ~ProbabilitySolver() { }

      /// @brief Computes the probability of meeting the deadlines.
      ///
      /// This function computes the steady state probability of meeting a 
      /// set of deadlines. It is supposed to have access to the internal 
      /// information of the task descriptor it applies to. The result of 
      /// the call is to fill in the proabilistic deadline structure within 
      /// the associated task descriptors.
      virtual void solve();

      /// @brief Check whether the solver has been executed.
      ///
      /// This function returns the state of the solved flag.
      ///
      /// @return the solved flag.
      bool is_solved() const {
        return solved;
      }

      ///@brief Resets the solver.
      ///
      /// This function resets the solver to make it ready for a new execution.
      virtual void reset() {

        /* Set the corresponding flags as unsolved */
        solved = false;
        pre_process_done = false;
        post_process_done = false;

      }

  };

}

#endif
