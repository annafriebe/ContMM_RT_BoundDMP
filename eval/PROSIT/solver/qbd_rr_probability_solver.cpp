/*!
 * @file    qbd_rr_probability_solver.cpp
 * 
 * @brief   This class defines an implementation of general functions for 
 *          QBD-based solvers. The core is the method for the computation 
 *          of the matrices. In addition we have some auxiliary function 
 *          for managing matrices.
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
#include "qbd_rr_probability_solver.hpp"

using namespace PrositAux;
using namespace Eigen;
using namespace std;

namespace PrositCore {

  /// @brief Prepare the computation.
  void QBDResourceReservationProbabilitySolver::preProcess() {

    uint64_t cmax, cmin;
    uint64_t zmax, zmin;
    uint64_t granularity;
    uint64_t Q;
    int m;                  ///< Number of modes.
    int64_t blockSize;      ///< Size of the submatrices.
    int64_t M;              ///< Size of the transition matrix.
    int64_t R0 = 0;         ///< # of repeated rows.
    int64_t R1 = 0;         ///< # of transient rows.
    int64_t H1;             ///< Total # of columns.
    int64_t H2;             ///< Total # of rows.

    /* Check whether the pre-process was already performed */
    if (pre_process_done) {

      if (verbose_flag) {
        cerr << "The pre-process phase was already performed for the task " 
             << task_descriptor->getName() << endl;
      }

      return;

    }

    /* Obtain the distribution of the computation time */
    unique_ptr<pmf> tmp(new pmf(*(task_descriptor->getComputationTime())[0]));
    unique_ptr<pmf> c = move(tmp->resample(task_descriptor->getGranularity(), true));

    /* Obtain information about the computation time */
    cmax = c->getMax();
    cmin = c->getMin();
    granularity = task_descriptor->getGranularity();

    /* Obtain the distribution of the interarrival time */
    unique_ptr<pmf> u(new pmf(*task_descriptor->getInterarrivalTime()));

    /* Obtain information about the interarrival time */
    zmax = u->getMax();
    zmin = u->getMin();

    /* Present online information */
    if (verbose_flag) {
      cerr << "Maximum computation time:  " << cmax << endl;
      cerr << "Minimum computation time:  " << cmin << endl;
      cerr << "Maximum interarrival time: " << zmax << endl;
      cerr << "Minimum interarrival time: " << zmin << endl;
    }

    /* Obtain the budget */
    Q = task_descriptor->getBudget();

    /* Present online information */
    if (verbose_flag) {

      /* The budget is a multiple of the granularity */
      if (Q % granularity) {
        cerr << "WARNING: The granularity is not a submultiple of the budget "
             "for task " << task_descriptor->getName() << endl;
      }

    }

    /* Adjust the budget */
    Q = Q / granularity;

    /* Obtain the number of modes */
    m = task_descriptor->getNumModes();
    MatrixXd mode_trans_prob = task_descriptor->getTransitionMatrix();

    /* Compute the boundaries of the transition matrix */
    R0 = ((zmin * Q) - cmin + 1) * m;
    R1 = (((zmax - zmin) * Q) * m);

    /* Check whether there are no repeated rows */
    if (R0 < 0) {
      R0 = 0;
    }

    H1 = R1 + ((cmax - cmin + 1) * m);
    H2 = R0 + R1;

    /* Obtain the maximum value */
    blockSize = max(H1, H2);
    M = 2 * blockSize;

    /* Present online information */
    if (verbose_flag) {
      cerr << "Repeated rows:             " << R0 << endl;
      cerr << "Transient rows:            " << R1 << endl;
      cerr << "Repeated colums:           " << H1 << endl;
      cerr << "Submatrices size:          " << blockSize << endl;
      cerr << "Matrix size:               " << M << endl;
      cerr << "Computing the transition matrix" << endl;
    }

    /* Initialize the probability transition matrix */
    MatrixXd P = MatrixXd::Zero(M, M);

    /* Initialize the matrix components */
    MatrixXd repeatedBlock = MatrixXd::Zero(m, blockSize);
    MatrixXd transientBlock = MatrixXd::Zero(m, blockSize);
    
    uint64_t idxCol = 0, idxRow = R0, timeIdx = 0;
    double s = 0.0, p_z = 0.0, p_c = 0.0;

    /* Construct the repeated block */
    for (int x = 0; x < m; x++) {

      /* Initialize the index for the column */
      idxCol = 0;

      /* Iterate over the length of the PMFs */
      for (uint64_t y = cmin; y <= cmax; y++) {

        /* Iterate over the arriving states */
        for (int z = 0; z < m; z++) {

          unique_ptr<pmf> tmp1(new pmf(*(task_descriptor->getComputationTime())[z]));
          unique_ptr<pmf> t = move(tmp1->resample(task_descriptor->getGranularity(), true));

          /* Compute the probability */
          repeatedBlock(x, idxCol) = mode_trans_prob(x, z) * (t->get(y - cmin));

          /* Update the index for the column */
          idxCol++;

        }

      }

    }

    /* Construct the part of the matrix that remains equal */
    for (int k = 0; k < (R0 / m); k++) {

      /* Add the repeated block to the matrix */
      P.block(k * m, 0, m, blockSize) = repeatedBlock;

    }

	  /* Construct the transient part */
	  /* Iterate over the current state */
    for (int h = (R0 / m); h < (H2 / m); h++) {

      for (int g = 0; g < m; g++) {

        idxCol = 0;

        for (int hp = 0; hp < H1; hp++) {

          for (int gp = 0; gp < m; gp++) {

            unique_ptr<pmf> tmp1(new pmf(*(task_descriptor->getComputationTime())[gp]));
            unique_ptr<pmf> t = move(tmp1->resample(task_descriptor->getGranularity(), true));

            /* Initialize the sum */
            s = 0.0;

            /* Iterate over the probability of inter-arrival time */
            for (uint64_t z = zmin; z <= zmax; z++) {

              /* Obtain the probability of the inter-arrival time */
              p_z = u->get(z - zmin);

              /* Compute the index for the computation time */
              timeIdx = (int)cmin + hp - max(0, (int)cmin + h - ((int)z * (int)Q));

              /* Obtain the probability of the computation time */
              if ((timeIdx >= cmin) && (timeIdx <= cmax)) {
                p_c = t->get(timeIdx - cmin);
              }
              else {
                p_c = 0;
              }

              /* Accumulate the probability */
              s += (p_z * mode_trans_prob(g, gp) * p_c);

            }

            P(idxRow, idxCol) = s;

            idxCol++;

          }

        }

        idxRow++;

      }
  
    }

    /* Select the shifting row */
    transientBlock = P.block(H2 - m, 0, m, blockSize);

    /* Initialize the column index */
    idxCol = m;

    /* The number of rows is greater than the number of columns */
    if (H2 >= H1) {

      /* Construct the shifting rows without cutting the pmf resampled vector */
      for (int k = (H2 / m); k < (M / m ); k++) {

        /* Add the repeated shifting block to the matrix */
        P.block(k * m, idxCol, m, blockSize) = transientBlock;

        /* Update the column index */
        idxCol = idxCol + m;

      }

    }
    else {

      /* Construct the shifting rows without cutting the pmf resampled vector */
      for (int k = (H2 / m); k < ((H2 / m) + ((M - blockSize) / m)); k++) {

        /* Add the repeated shifting block to the matrix */
        P.block(k * m, idxCol, m, blockSize) = transientBlock;

        /* Update the column index */
        idxCol = idxCol + m;

      }

      /* Number of elements to be removed from the pmf resampled vector */
      int n = 1;

      /* Construct the rows that shift cutting the pmf resampled vector */
      for (int k = ((H2 / m) + ((M - blockSize) / m)); k < (M / m); k++) {

        /* Add the cut shifting block to the matrix */
        P.block(k * m, idxCol, m, blockSize - (n * m)) = transientBlock.block(0, 0, m, 
            blockSize - (n * m));

        /* Update the column index */
        idxCol = idxCol + m;

        /* Update the number of columns to cut*/
        n++;

      }

    }

    /* Present online information */
    if (verbose_flag) {
      cerr << "Extracting submatrices" << endl;
    }

    /* Initialize submatrices */
    A0 = MatrixXd(blockSize, blockSize);
    B0 = MatrixXd(blockSize, blockSize);
    A1 = MatrixXd(blockSize, blockSize);
    A2 = MatrixXd(blockSize, blockSize);

    /* Extract submatrices */
    extractSubmatrices(P, blockSize, B0, A0, A1, A2);

    /* Initialize matrix */
    R = MatrixXd(A0.rows(), A0.cols());

  }

  /// @brief Post processing after the solution.
  void QBDResourceReservationProbabilitySolver::postProcess() {

    /* Compute the stationary vector */
    if (!computePi0()) {

      /* Present online information */
      if (task_descriptor->getVerboseFlag()) {
        cerr << "WARNING: There were anomalies in the computation of pi0" 
             << endl;
      }

    }

  }

  /// @brief Fills in the probability map.
  void QBDResourceReservationProbabilitySolver::fillInProbabilityMap() {

    //cout << pi0.transpose() << endl;
    uint64_t cmin;
    uint64_t Ts, Q;
    int m;
    double prob = 0.0;

    /* Obtain the distribution of the computation time */
    unique_ptr<pmf> tmp(new pmf(*(task_descriptor->getComputationTime())[0]));
    unique_ptr<pmf> c = move(tmp->resample(task_descriptor->getGranularity(), true));

    /* Obtain information about the computation time */
    cmin = c->getMin();

    /* Obtain the scheduling parameters */
    Ts = task_descriptor->getServerPeriod();
    Q = task_descriptor->getBudget() / task_descriptor->getGranularity();
    m = task_descriptor->getNumModes();

    /* Obtain the deadline probability map */
    DeadlineProbabilityMap *pm = task_descriptor->getProbabilisticDeadlines();
    DeadlineProbabilityMapIter pmi;

    /* Obtain the maximum number of iterations */
    int H = ceil((double(pm->size()) * double(Q)) / (double(pi0.size()) / double(m)));

    RowVectorXd pi = pi0;
    
    /* Iterate for the maximum length of the vector */
    for (int h = 0; h < H; h++) {

      /* Iterate over the stationary vector */
      for (int i = 0; i < (pi.size() / m); i++) {

        for (int k = 0; k < m; k++) {

          /* Accumulate the probability */
          prob += pi((i * m) + k);

        }

        /* The iteration has reached a probabilistic deadline */
        if (((cmin + i + (h * (pi.size() / m))) % Q) == 0) {

          uint64_t delta = ((cmin + i + (h * (pi.size() / m))) / Q) * Ts;

          /* Obtain the probabilisic deadline */
          if ((pmi = pm->find(delta)) != pm->end()) {

            /* Associate the probability to the probabilistic deadline */
            if (prob <= 1.0) {
              (*pmi).second = prob;
            }
            else {
              (*pmi).second = 1.0;
            }

          }

        }

      }

      /* Obtain the next chunk of probabilities */
      pi = pi * R;

    }

  }

  /// @brief Computes the stationary vector pi0.
  bool QBDResourceReservationProbabilitySolver::computePi0() {

    /* Check whether there is a linked task */
    if (!task_descriptor) {
      EXC_PRINT("ERROR: The computation of pi0 requires a task descriptor.");
    }

    /* Check whether thw solver has been applied */
    if (!solved) {
      EXC_PRINT_2("ERROR: The solver has not been applied for task ", 
          task_descriptor->getName());
    }

    /* Check the size of the matrix R */
    if (R.rows() != R.cols()) {
      EXC_PRINT("ERROR: The matrix R has to be square");
    }

    /* The stationary vector was already computed */
    if (post_process_done && verbose_flag) {

      cerr << "WARNING: Computation of pi0 required twice. Task " 
           << task_descriptor->getName() << endl;

    }

    /* Check the proper dimension of the matrices */
    if (!check_sizes(A0, R) || !check_sizes(A0, B0)) {
      EXC_PRINT("ERROR: The matrices A0, A1 and A2 have to be square and "
          "equal size");
    }

    /* Negative coefficients in matrix R */
    if ((R.minCoeff() < 0) && verbose_flag) {
      cerr << "WARNING: The matrix R has negative coeeficients" << endl;
    }

    /* Compute the eigenvalues and eigenvectors of matrix R */
    EigenSolver<MatrixXd> eigensolver(R);

    /* The eigenvalues were correctly computed */
    if (eigensolver.info() != Success) {

      /* Present online information */
      if (verbose_flag) {
        cerr << "WARNING: It was not possible to compute the eigenvalues of R" 
             << endl;
      }

      return false;

    }

    /* The spectral radius is greater than 1 */
    if (((eigensolver.eigenvalues()).array().abs().maxCoeff() > 1) && verbose_flag) {
      EXC_PRINT("ERROR: The matrix R has spectral radius greater than 1");
    }

    int n = R.rows();
    MatrixXd Id = MatrixXd::Identity(n, n);
    VectorXd u = VectorXd::Ones(n);
    MatrixXd M(n, n + 1);
    M.block(0, 0, n, n) = B0 + R * A0 - Id;

    /* ******************************************************************* */
    /* MATLAB algorithm */
    MatrixXd temp = (Id - R).inverse();
    RowVectorXd pi = stat(B0 + R * A0);
    pi0 = pi / (pi * temp * u);

/*    if (pi0.minCoeff() < 0.0) {
      EXC_PRINT("ERROR: The stationary vector pi0 has negative elements");
    }*/

    /* Compute the complete vector
    double sumpi = pi.sum();
    long numit = 0;

    while ((sumpi < (1-10^(-10))) && (numit < 500)) {

      pi(numit+1,1:m)=pi(numit,:)*R; % compute pi_(numit+1)
      numit++;
      sumpi=sumpi+sum(pi(numit,:));

    }
    */
    /* ******************************************************************* */

    /* ******************************************************************* */
    /* PROSIT algorithm: Extremely Time Consuming
    M.block(0, n, n, 1) = (Id - R).inverse() * u;
    FullPivLU<MatrixXd> lu_decomp(M);

    if (lu_decomp.rank() < n) {

      if (verbose_flag)
        cerr << "WARNING: There is no unique solution" << endl;

      return false;

    }

    RowVectorXd work(n + 1);
    work.setZero();
    work(n) = 1;
    MatrixXd W1;
    pseudoInverse<MatrixXd>(M, W1);

    pi0 = work * W1;*/

    /* ******************************************************************* */

    return true;

  }

  /// @brief Extract the QBDP submatrices.
  void QBDResourceReservationProbabilitySolver::extractSubmatrices(const MatrixXd &P,
      uint64_t dim, MatrixXd &B, MatrixXd &A0, MatrixXd &A1, MatrixXd &A2) {

    B = P.block(0, 0, dim, dim);
    A0 = P.block(dim, 0, dim, dim);
    A1 = P.block(dim, dim, dim, dim);
    A2 = P.block(0, dim, dim, dim);

  }

}
