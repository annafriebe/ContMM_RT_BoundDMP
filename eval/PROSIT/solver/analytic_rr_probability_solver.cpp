/*!
 * @file    analytic_rr_probability_solver.cpp
 * 
 * @brief   This class defines an implementation for the Analytic solver for 
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
#include "analytic_rr_probability_solver.hpp"

using namespace PrositAux;
using namespace Eigen;
using namespace std;

namespace PrositCore {

  void AnalyticResourceReservationProbabilitySolver::computeNextVector(RowVectorXd &pi) {

    uint64_t size = pi.cols();

    if (size >= 2) {
      RowVectorXd tmp_pi = pi.block(0, 1, 1, size - 1);
      VectorXd tmp_alfa(size - 1);
      tmp_alfa << (-1 * Alfa.reverse().transpose()), (Alfa.sum() + 1.0);
    
      for (uint64_t i = 0; i < size; i++) {

        pi(i) = tmp_pi * tmp_alfa;
        tmp_pi << tmp_pi.block(0, 1, 1, size - 2), pi(i);

      }

    }
    
  }

  /// @brief Applies the algorithm.
  void AnalyticResourceReservationProbabilitySolver::applyAlgorithm() {

    double a0 = 0.0;
    double pi_0 = 1.0;

    /* Obtain the PMF of the computation time */
    unique_ptr<pmf> tmp(new pmf(*(task_descriptor->getComputationTime())[0]));
    unique_ptr<pmf> p = move(tmp->resample(task_descriptor->getGranularity(), true));

    /* Obtain information about the computation time */
    int64_t cmin = int64_t(p->getMin());
    int64_t cmax = int64_t(p->getMax());

    /* Obtain the assigned CPU time */
    int64_t bw = (task_descriptor->bandwidth * task_descriptor->getPeriod()) / task_descriptor->getGranularity();

    /* Obtain the size of the pi0 vector */
    int64_t size = int64_t(p->getSize()) - bw + cmin + 1;

    /* Obtain the denominator of the analytic bound */
    for (int64_t i = 0; i < (bw - cmin); i++) {
      a0 += p->get(i);
    }

    /* Present online information */
    if (verbose_flag) {
      cerr << "Bandwidth: " << bw << endl;
      cerr << "WCET: " << cmax << endl;
      cerr << "Preparing analytic form" << endl;
    }

    /* The assigned computation time is enough for the task */
    if (bw > cmax) {

      /* Present online information */
      if (verbose_flag) {
        cerr << "The bandwidth is greater than the worst case requirements" << endl;
      }

      pi_0 = 1.0;

    }

    if (abs(a0) < 1e-10) {

      /* Present online information */
      if (verbose_flag) {
        cerr << "The bandwidth is too small" << endl;
      }

      pi_0 = 0.0;

    }
    else {

      /* Iterate over the valid computation times */
      for (int64_t i = (bw + 1); i <= cmax; i++) {

        /* Update the probability of respecting the deadline */
        pi_0 -= (i - bw) * (p->get(i - cmin) / a0);

      }

    }

    /* Present online information */
    if (verbose_flag) {
      cerr << "The analytic computation is completed" << endl;
    }

    /* Check whether the probability is smaller than zero */
    if (pi_0 < 0) {
      pi_0 = 0;
    }

    /* Build the complete probability vector */
    pi0 = RowVectorXd(size);
    pi0(0) = pi_0;

    /* Compute the row vector of alpha's */
    if (size >= 2) {

      Alfa = RowVectorXd(size - 2);
      for (int64_t i = 0; i < (size - 2); i++) {

        Alfa(i) = p->get(bw - cmin + i + 1) / a0;

      }

      /* Build the Equation 27 for pi_1 */
      pi0(1) = Alfa.sum() * pi_0;
    
      double term = Alfa.sum() + 1.0;

      for (int64_t l = 2; l < size; l++) {

        double pos_term = term * pi0(l - 1);

        double neg_term = 0.0;

        int64_t idx = min(size - 1, l);

        for (int64_t j = 2; j <= idx; j++) {

          neg_term += Alfa(j - 2) * pi0(idx - j);

        }

        pi0(l) = pos_term - neg_term;

      }
    }

    return;

  }

  /// @brief Fills in the probability map.
  void AnalyticResourceReservationProbabilitySolver::fillInProbabilityMap() {

    uint64_t cmin;
    uint64_t Ts, Q, delta;
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
    delta = task_descriptor->getPeriod();

    /* Obtain the deadline probability map */
    DeadlineProbabilityMap *pm = task_descriptor->getProbabilisticDeadlines();
    DeadlineProbabilityMapIter pmi;

    RowVectorXd pi = pi0;

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

          /* Update the deadline */
          delta += Ts;

        }

      }

      /* Obtain the next chunk of probabilities */
      computeNextVector(pi);

    }

    return;

  }

}
