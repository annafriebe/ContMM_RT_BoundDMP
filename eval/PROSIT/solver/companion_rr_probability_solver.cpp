/*!
 * @file    companion_rr_probability_solver.hpp
 * 
 * @brief   This class defines an implementation for the Companion form solver
 *          for resource reservation tasks described as a Quasi-Birth-Death 
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
#include "companion_rr_probability_solver.hpp"

using namespace PrositAux;
using namespace Eigen;
using namespace std;

namespace PrositCore {

  complex<double> CompanionResourceReservationProbabilitySolver::GammaFunctionsTotal(const complex<double> &b, uint64_t dim, const RowVectorXd &Alpha) {

    size_t size = Alpha.size();
    complex<double> gamma = pow<double>(b, dim);

    for (uint64_t i = 1; i <= dim; i++) {
      gamma += pow<double>(b, int(dim) - int(i)) * -Alpha(size - i);
    }

    return gamma;

  }

  /// @brief Check the conditions to start.
  bool CompanionResourceReservationProbabilitySolver::checkList() {

    /* Check the general conditions */
    bool ok = ResourceReservationProbabilitySolver::checkList();

    /* The solver was already applied */
    if (!ok) {

      return false;

    }

    /* Check the periodicity of the task */
    if (!(task_descriptor->isPeriodic())) {

      if (verbose_flag) {
        cerr << "The companion solver currently applies to periodic single mode tasks." << endl;
      }

      return false;

    }

    return true;

  }

  /// @brief Prepare the computation.
  void CompanionResourceReservationProbabilitySolver::preProcess() {

    /* Check whether the solver has been linked to a task */
    if (!linked_flag) {

      EXC_PRINT("ERROR: The solver has not been linked to a task.");

    }

    /* Check whether the solver was already called */
    if (solved) {

      if (verbose_flag)
        cerr << "The solver was already called for the task "
             << task_descriptor->getName() << endl;

      return;

    }

    /* Check whether the pre-process was already performed */
    if (pre_process_done) {

      if (verbose_flag) {
        cerr << "The pre-process phase was already performed for the task " 
             << task_descriptor->getName() << endl;
      }

      return;

    }

    /* Obtain the PMF of the computation time */
    unique_ptr<pmf> tmp(new pmf(*(task_descriptor->getComputationTime())[0]));
    unique_ptr<pmf> p = move(tmp->resample(task_descriptor->getGranularity(), true));

    /* Obtain the sheduling parameters */
    uint64_t Q = task_descriptor->getBudget() / task_descriptor->getGranularity();
    uint64_t N = task_descriptor->getPeriod() / task_descriptor->getServerPeriod();

    /* Obtain information of the PMF */
    uint64_t H = (N * Q) - p->getMin();

    uint64_t size = p->getSize() - 1;

    /* Obtain the value of a0 */
    double a0 = p->get(0);

    /* Compute the row vector of alpha's */
    Alfa = RowVectorXd(size);
    for (uint64_t i = 1; i <= size; i++) {

      double ai = p->get(i);

      if (i == H) {
        Alfa(size - i) = (1 - ai) / a0;
      }
      else {
        Alfa(size - i) = - ai / a0;
      }

    }

    /* Define the auxiliary elements of the companion matrix */
    VectorXd z = VectorXd::Zero(size - 1);
    MatrixXd I = MatrixXd::Identity(size - 1, size - 1);

    /* Build the companion matrix */
    W = MatrixXd(size, size);
    W.block(0, 0, size - 1, 1) = z;
    W.block(0, 1, size - 1, size - 1) = I;
    W.block(size - 1, 0, 1, size) = Alfa;

  }

  /// @brief Executes the iterations of the Companion form method.
  void CompanionResourceReservationProbabilitySolver::applyAlgorithm() {

    double pi_0;
    uint64_t size = W.rows();

    /* Obtain the eigenvalues of the companion matrix */
    EigenSolver<MatrixXd> eigensolver(W.transpose());

    /* Check the correctness of the eigensolver */
    if (eigensolver.info() != Success) {
      EXC_PRINT("ERROR: Cannot compute the eigenvalues of the companion matrix.");
    }

    /* Obtain the elements to satisfy the equation A^n = P * B^n * P^{-1}*/
    MatrixXcd P, Pinv;
    P = eigensolver.eigenvectors();
    pseudoInverse<MatrixXcd>(P, Pinv);

    /* Compute the power matrix (W^j) */
    D = (eigensolver.eigenvalues().array().pow(size)).matrix().asDiagonal();

    /* Variable for the computation of pi0 */
    complex<double> product(1, 0);

    /* Variable for the indices of the unstable eigenvectors */
    vector<int> unst_idx;

    /* Iterate over the eigenvalues */
    for (uint64_t i = 0; i < size; i++) {

      /* Obtain the corresponding eigenvalue */
      complex<double> p_i = eigensolver.eigenvalues()[i];

      /* The stable eigenvalues are used for the computation of pi0 */
      if (abs(p_i) < (1 - 1e-10)) {

        product *= (complex<double>(1, 0) - p_i);

      }
      /* The unstable eigenvalues are stored for other computations */
      else {

        unst_idx.push_back(i);
        num_unst_eig++;
        D(i, i) = 0;

      }

    }

    /* Obtain the power matrix W^j */
    D = (P * D * Pinv).transpose();

    /* Check whether there is an imaginary part in the probability */
    if (product.imag() > 1e-14) {
      EXC_PRINT("ERROR: The probability did not result into a real number.");
    }

    /* Obtain the value of pi0 */
    if ((product.real() < 0.0) || (product.real() > (1 + 1e-14))) {
      EXC_PRINT("ERROR: The companion solver did not produce the correct result.");
    }
    else {
      pi_0 = product.real();
    }

    /* Variables for the unstable eigenvectors (eigenvalue >= 1) */
    MatrixXcd A = MatrixXcd::Zero(size - 1, size - 1);
    VectorXcd b = VectorXcd::Zero(size - 1);

    /* Iterate over the set of unstable eigenvalues */
    for (uint64_t h = 0; h < num_unst_eig; h++) {

      /* Obtain the eigenvectors in a row fashion */
      RowVectorXcd Ah = eigensolver.eigenvectors().col(unst_idx[h]).transpose();

      /* Build the system of equations (Ax = b) */
      A.block(h, 0, 1, size - 1) = Ah.block(0, 1, 1, size - 1);
      b(h) = -Ah(0) * pi_0;

    }

    /* Build the Equation 5 for l = 0 */
    for (uint64_t h = 0; h < num_unst_eig - 1; h++) {

      A(num_unst_eig, num_unst_eig - 2 - h) = -GammaFunctionsTotal(complex<double>(1, 0), h + 1, Alfa);

    }

    /* Define the coefficient for pi_H */
    A(num_unst_eig, num_unst_eig - 1) = -1;

    /* Obtain the coefficient without unknown variable */
    for (uint64_t h = num_unst_eig + 1; h <= size; h++) {
      b(num_unst_eig) += Alfa(size - h) * pi_0;
    }

    /* Obtain the number of missing equations */
    uint64_t num_eq = size - num_unst_eig - 2;

    /* Build the Equation 5 for l >= 1 */
    for (uint64_t l = 0; l < num_eq; l++) {

      /* Define the coefficients for the unknown variables */
      A.block(num_unst_eig + l + 1, 0, 1, num_unst_eig + l) = Alfa.block(0, size - num_unst_eig - l, 1, num_unst_eig + l).cast<complex<double>>();

      /* Define the coefficient for pi_{H + l + 1} */
      A(num_unst_eig + l + 1, num_unst_eig + l) = complex<double>(-1, 0);

      /* Obtain the coefficient without unknown variable */
      b(num_unst_eig + l + 1) = -Alfa(size - num_unst_eig - l - 1) * pi_0;

    }

    /* Solve the system of equations */
    MatrixXcd Ainv;
    pseudoInverse<MatrixXcd>(A, Ainv);
    VectorXcd x = Ainv * b;

    /* Build the complete probability vector */
    pi0 = RowVectorXd(size);
    pi0 << pi_0, x.real().transpose();

    return;

  }

  /// @brief Fills in the probability map.
  void CompanionResourceReservationProbabilitySolver::fillInProbabilityMap() {

    uint64_t cmin;
    uint64_t Ts, Q, idx_cut = 0;
    double prob = 0.0;
    uint64_t size = pi0.cols();

    /* Obtain the distribution of the computation time */
    unique_ptr<pmf> tmp(new pmf(*(task_descriptor->getComputationTime())[0]));
    unique_ptr<pmf> c = move(tmp->resample(task_descriptor->getGranularity(), true));

    /* Obtain information about the computation time */
    cmin = c->getMin();

    /* Obtain the scheduling parameters */
    Ts = task_descriptor->getServerPeriod();
    Q = task_descriptor->getBudget() / task_descriptor->getGranularity();

    /* Obtain the deadline probability map */
    DeadlineProbabilityMap *pm = task_descriptor->getProbabilisticDeadlines();
    DeadlineProbabilityMapIter pmi;

    RowVectorXd pi = RowVectorXd::Zero(size);

    /* Compute the probabilities below the task period */
    for (uint64_t i = 0; i < size; i++) {

      /* Obtain the value for the probability vector */
      for (uint64_t k = 0; k <= i; k++) {

        pi(i) += c->get(i - k) * pi0(k);

      }

      /* Accumulate the probability */
      prob += pi(i);

      /* The vector has reached the task period */
      if (abs(prob - pi0(0)) < 1e-14) {
        idx_cut = i;
      }

    }

    prob = 0.0;

    /* Obtain the maximum number of iterations */
    int H = ceil((double(pm->size()) * double(Q)) / double(size));

    /* Iterate for the maximum length of the vector */
    for (int h = 0; h < H; h++) {

      /* Iterate over the stationary vector */
      for (uint64_t i = 0; i < size; i++) {

        /* Accumulate the probability */
        prob += pi(i);

        /* The iteration has reached a probabilistic deadline */
        if (((cmin + i + (h * size)) % Q) == 0) {

          uint64_t delta = ((cmin + i + (h * size)) / Q) * Ts;

          /* Obtain the probabilistic deadline */
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

      /* Obtain the next chunk of probabilities pi^j = W^j * pi^0 */
      RowVectorXd tmp_pi = ((D * pi0.transpose()).transpose()).real();

      /* Build the vector for the computation of the probability */
      pi << pi0.block(0, size - idx_cut, 1, idx_cut), tmp_pi.block(0, 0, 1, size - idx_cut);

      /* Update the probability vector for the next iteration */
      pi0 = tmp_pi;

    }

    return;

  }

}
