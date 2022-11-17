/*!
 * @file    cr_qbd_rr_probability_solver.cpp
 * 
 * @brief   This class defines the implementation of the Cyclic Reduction [1]
 *          solver for resource reservation tasks based on algorithms for 
 *          numeric solution of Quasi-Birth-Death processes.
 * 
 *          [1] B. Meing and D. Bini Improved cyclic reduction for solving 
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
#include "cr_qbd_rr_probability_solver.hpp"

using namespace PrositAux;
using namespace Eigen;
using namespace std;

namespace PrositCore {

  /// @brief Executes the iterations of the Cyclic Reduction method.
  void CRQBDResourceReservationProbabilitySolver::applyAlgorithm() {

    /* Check the size of the matrices */
    if (!check_sizes(A0, A1) || !check_sizes(A1, A2) || !check_sizes(A1, R)) {
      EXC_PRINT("ERROR: The matrices A0, A1 and A2 have to be square and "
          "equal size");
    }

    /* Initialize the auxiliary matrices */
    MatrixXd A0c = A0;
    MatrixXd A1c = A1;
    MatrixXd A2c = A2;
    double drift = 0;

    MatrixXd Id = MatrixXd::Identity(A2c.rows(), A2c.cols());
    VectorXd D = A2c.diagonal();
    VectorXd u = VectorXd::Ones(A1c.rows());
    RowVectorXd uT = RowVectorXd::Ones(A1c.cols());
    uT = 1 / (A2c.cols()) * uT;

    /* Check the use of discrete time systems */
    if (D.sum() < 0) {
      EXC_PRINT("ERROR: Only discrete time supported at the moment");
    }

    if (shift_flag) {

      RowVectorXd theta = stat(A2c + A1c + A0c);
      double tmp;

      drift = theta * A2c.rowwise().sum();
      tmp = theta * A0c.rowwise().sum();
      drift = drift - tmp;

      if (drift < 0) { // MC is transient: use its dual
        A0c = A0c - u * (theta * A0c);
        A1c = A1c + u * (theta * A2c);
      }
      else {
        A2c = A2c - A2c.rowwise().sum() * uT;
        A1c = A1c + A2.rowwise().sum() * uT;
      }

    }

    /* Start the cyclic reduction algorithm */
    MatrixXd A = A1c;
    MatrixXd B = A2c;
    MatrixXd C = A0c;
    MatrixXd Ahat = A;
    
    double check = 1.0;
    long numit = 0;

    while ((check > 1e-14) && numit < max_iter) {

      MatrixXd Atemp;
      MatrixXd BAtemp;

      Atemp = (Id - A).inverse();
      BAtemp = B * Atemp;
      Atemp = C * Atemp;
      Ahat = Ahat + BAtemp * C;
      A = A + BAtemp * C + Atemp * B;
      B = BAtemp * B;
      C = Atemp * C;
      numit++;

      check = min(InfinityNorm<MatrixXd>(B), InfinityNorm<MatrixXd>(C));

      /* Present online information */
      if (verbose_flag) {
        cerr << "After " << numit << " iterations " << check << " reached" << endl;
      }

    }

    /* The maximum number of iterations has been reached */
    if (numit == max_iter) {
      cerr << "WARNING: The maximum number of iterations has been reached " << endl;
    }

    /* Compute the G matrix */
    G = ((Id - Ahat).inverse()) * A0c;

    if (shift_flag) {

      if (drift < 0) {
        A1c = A1;
        A0c = A2;
      }
      else {
        G = G + u * uT;
        A1c = A1;
        A2c = A0;
      }

    }

    /* Present online information */
    if (verbose_flag) {
      cerr << "Final Residual Error for G: " 
           << InfinityNorm<MatrixXd>(G - A0c - ((A1c + (A2c * G)) * G)) << endl;
    }

    /* Compute the R matrix */
    R = A2c * (Id - (A1c + A2c * G)).inverse();

    /* Present online information */
    if (verbose_flag) {
      cerr << "Final Residual Error for R: " 
           << PrositAux::InfinityNorm<MatrixXd>(R - A2c - R * (A1c + R * A0c)) 
           << endl;
    }

    /* Compute the U matrix */
    U = A1c + R * A0c;

    /* Present online information */
    if (verbose_flag) {
      cerr << "Final Residual Error for U: " 
           << PrositAux::InfinityNorm<MatrixXd>(U - A1c - A2c * (Id - U).inverse() * A0c) 
           << endl;
    }

  }

}
