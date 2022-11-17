/*!
 * @file    probability_solver.cpp
 * 
 * @brief   This class contains the general implementation for probability 
 *          solvers.
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
#include "probability_solver.hpp"

using namespace std;

namespace PrositCore {

  /// @brief Computes the probability of meeting the deadlines.
  void ProbabilitySolver::solve() {

    /* Check the conditions to start */
    if (!checkList()) {

      /* There were error in the initial conditions */
      if (linked_flag && verbose_flag) {
        cerr << "Probability solver will not execute." << endl;
      }

      return;

    }

    /* Perform the pre-process phase */
    preProcess();

    /* The pre-processing was successfully completed */
    pre_process_done = true;
    if (verbose_flag) {
      cerr << "Pre-processing completed." << endl;
    }

    /* Apply the probability solver */
    applyAlgorithm();

    /* The algorithm was successfully applied */
    solved = true;
    if (verbose_flag) {
      cerr << "Solver computation completed" << endl;
    }

    /* Perform the post-process phase */
    postProcess();

    /* The post-processing was successfully completed */
    post_process_done = true;
    if (verbose_flag) {
      cerr << "Post-processing completed" << endl;
    }

    /* Compute the probabilistic deadlines */
    fillInProbabilityMap();

    /* The probability map was successfully filled in */
    if (verbose_flag) {
      cerr << "Probability map correctly filled" << endl;
    }

  }

}
