/*!
 * @file    lt_qbd_rr_probability_solver.cpp
 * 
 * @brief   This class defines the implementation of the Latouche Reduction [1]
 *          solver for resource reservation tasks based on algorithms for 
 *          numeric solution of Quasi-Birth-Death processes.
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
#include "lt_qbd_rr_probability_solver.hpp"

using namespace PrositAux;
using namespace Eigen;
using namespace std;

namespace PrositCore {

  /// @brief Executes the iterations of the Latouche method.
  void LatoucheQBDResourceReservationProbabilitySolver::applyAlgorithm() {

    MatrixXd Rnew = A2;
    R.setZero();
    int cnt = 0;

    /* Check the size of the matrices */
    if (!check_sizes(A0, A1) || !check_sizes(A1, A2) || !check_sizes(A1, R)) {
      EXC_PRINT("ERROR: The matrices A0, A1 and A2 have to be square and "
          "equal size");
    }

    /* Repeat the algorithm the required number of times */
    while ((InfinityNorm(R - Rnew) > epsilon) && (cnt < max_iter)) {

      /* Compute the solution for R */
      R = Rnew;
      Rnew = A2 + R * A1 + R * R * A0;

      /* Update the number of iterations */
      cnt++;

    }

    /* The maximum number of iterations has been reached */
    if (cnt == max_iter) {
      cerr << "WARNING: The maximum number of iterations has been reached." << endl;
    }

    /* Present online information */
    if (verbose_flag) {
      cerr << "Final Error for R: " << InfinityNorm(R - Rnew) << endl;
    }

    /* Update the final value for the R matrix */
    R = Rnew;

  }

}
