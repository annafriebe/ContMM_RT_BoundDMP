/*!
 * @file    fp_probability_solver.cpp
 * 
 * @brief   This class contains the implementation for generic probability
 *          solvers for a set of fixed priority tasks. In essence, a
 *          probability solver wraps the different algorithms used to compute
 *          probabilistic deadlines. The solver for fixed priority tasks
 *          considers sets of tasks.
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
#include "fp_probability_solver.hpp"

using namespace PrositAux;
using namespace Eigen;
using namespace std;

namespace PrositCore {

  void FixedPriorityProbabilitySolver::shiftDistribution(unique_ptr<pmf> &dist, uint64_t delta) {

    if ((delta != 0) && ((int64_t(dist->getMin()) - int64_t(delta)) >= 0)) {

      dist->setMin(dist->getMin() - delta);
      dist->setMax(dist->getMax() - delta);

    }

  }

  void FixedPriorityProbabilitySolver::normalizeDistribution(unique_ptr<pmf> &dist) {

    /* Obtain the current sum of the distribution */
    double summation = dist->sum();

    /* Iterate over the distribution */
    for (uint64_t k = dist->getMin(); k <= dist->getMax(); k++) {

      /* Normalize each element of the distribution */
      dist->set(k - dist->getMin(), dist->get(k - dist->getMin()) / summation);

    }

  }

  unique_ptr<pmf> FixedPriorityProbabilitySolver::addDistribution(unique_ptr<pmf> &first, const unique_ptr<pmf> &second) {

    /* Obtain the limits of the summation */
    uint64_t kmin = (first->getMin() <= second->getMin()) ? first->getMin() : second->getMin();
    uint64_t kmax = (first->getMax() >= second->getMax()) ? first->getMax() : second->getMax();

    /* Variable for the summation of the distributions */
    unique_ptr<pmf> summation(new pmf(kmin, kmax, 1e-4));

    /* Iterate over the elements of the final sum */
    for (uint64_t k = summation->getMin(); k <= summation->getMax(); k++) {

      double summ = 0.0;

      /* Add the corresponding elements of the first distribution */
      if ((k >= first->getMin()) && (k <= first->getMax())) {
        summ += first->get(k - first->getMin());
      }
      
      /* Add the corresponding elements of the second distribution */
      if ((k >= second->getMin()) && (k <= second->getMax())) {
        summ += second->get(k - second->getMin());
      }

      /* Set the final summation */
      summation->set(k - summation->getMin(), summ);

    }

    return summation;

  }

  unique_ptr<pmf> FixedPriorityProbabilitySolver::splitConvolve(unique_ptr<pmf> &first, const unique_ptr<pmf> &second, uint64_t delta_second) {

    /* Size of the distributions starting from 0 */
    uint64_t size1 = first->getMax() + 1;

    if (size1 < delta_second) {

      unique_ptr<pmf> clone(new pmf(first->getMin(), first->getMax(), 1e-4));

      for (uint64_t k = 0; k < first->getSize(); k++) {

        clone->set(k, first->get(k));

      }

      return clone;

    }

    /* Variable for the tail of the distribution */
    unique_ptr<pmf> tail(new pmf(delta_second + 1, first->getMax(), 1e-4));

    /* Split the original distribution to obtain the tail */
    for (uint64_t k = delta_second + 1; k <= first->getMax(); k++) {

      tail->set(k - (delta_second + 1), first->get(k - first->getMin()));

    }

    /* Perform the convolution between the tail and the computation time */
    unique_ptr<pmf> convolved = move(convolutionsConvolve(tail, second));

    /* Variable for the merged distribution */
    unique_ptr<pmf> merged(new pmf(first->getMin(), convolved->getMax(), 1e-4));

    /* Merge the head part of the distribution */
    for (uint64_t k = 0; k <= (delta_second - first->getMin()); k++) {

      merged->set(k, first->get(k));

    }

    /* Merge the convolved part of the distribution */
    for (uint64_t k = convolved->getMin(); k <= convolved->getMax(); k++) {

      merged->set(k - merged->getMin(), convolved->get(k - convolved->getMin()));

    }

    return merged;

  }

  /// @brief Computes the response time for a job.
  unique_ptr<pmf> FixedPriorityProbabilitySolver::getJobResponseTime(unique_ptr<pmf> &job_back, uint64_t id_task, uint64_t curr_activation) {

    /* Obtain the schedule of activations */
    vector<pair<uint64_t, uint64_t>> sched = task_schedule->getSchedule();

    /* Obtain the release time and the deadline of the current activation */
    uint64_t release_time = sched[curr_activation].second;
    uint64_t abs_deadline = ((task_schedule->getTaskSet())[id_task])->getDeadlineStep() + release_time;

    unique_ptr<pmf> resp_time = move(convolveShrink(job_back, ((task_schedule->getTaskSet())[id_task])->getComputationTime(), 0));

    /* Iterate over the upcoming activations */
    for (size_t i = curr_activation + 1; (sched[i].second < abs_deadline) && i < sched.size(); i++) {

      if (task_schedule->hasPriority(sched[i].first, id_task)) {

        resp_time = splitConvolve(resp_time, ((task_schedule->getTaskSet())[sched[i].first])->getComputationTime(), sched[i].second - release_time);

      }

    }

    return resp_time;

  }

  /// @brief Perform the convolution of two given distributions.
  unique_ptr<pmf> FixedPriorityProbabilitySolver::convolutionsConvolve(const unique_ptr<pmf> &first, const unique_ptr<pmf> &second) {

    /* Size of the distributions starting from 0 */
    uint64_t size1 = first->getMax() + 1;
    uint64_t size2 = second->getMax() + 1;

    /* Variables for the distributions ranges */
    uint64_t cmin = numeric_limits<uint64_t>::max();
    uint64_t cmax = numeric_limits<uint64_t>::min();

    /* Temporal distribution to hold the convolution */
    unique_ptr<pmf> convolved(new pmf(0, size1 + size2 - 1, 1e-4));

    /* Iterate over the possible values of the convolution */
    for (uint64_t n = 0; n < (size1 + size2 - 1); n++) {

      uint64_t kmin, kmax, k;
      double acum = 0.0;

      /* Obtain the limits of the convolution */
      kmin = (n >= size2 - 1) ? n - (size2 - 1) : 0;
      kmax = (n < size1 - 1) ? n : size1 - 1;

      /* Iterate over the limits to perform the convolution */
      for (k = kmin; k <= kmax; k++) {

        double term1, term2;

        /* Verify whether the desired index belongs to the distribution */
        if ((k >= first->getMin()) && (k <= first->getMax())) {
          term1 = first->get(k - first->getMin());
        }
        else {
          term1 = 0.0;
        }

        /* Verify whether the desired index belongs to the distribution */
        if (((n - k) >= second->getMin()) && ((n - k) <= second->getMax())) {
          term2 = second->get(n - k - second->getMin());
        }
        else {
          term2 = 0.0;
        }

        /* Accumulate the probability */
        acum += (term1 * term2);

      }

      /* Set the values of the distribution */
      convolved->set(n, acum);

      /* Obtain the proper range of the convolved distribution */
      if (acum != 0) {

        /* Determine the minimum index */
        if (n < cmin) {
          cmin = n;
        }

        /* Determine the maximum index */
        if (n > cmax) {
          cmax = n;
        }

      }

    }

    /* Proper convolved distribution */
    unique_ptr<pmf> clone(new pmf(cmin, cmax, 1e-4));

    /* Iterate in the range of the new distribution */
    for (uint64_t i = cmin; i <= cmax; i++) {

      /* Set the values of the distribution */
      clone->set(i - cmin, convolved->get(i));

    }

    return clone;

  }

  /// @brief Perform the shrinking of a given distributions
  unique_ptr<pmf> FixedPriorityProbabilitySolver::convolutionsShrink(unique_ptr<pmf> &dist, uint64_t delta) {

    /* Size of the distribution starting from 0 */
    uint64_t size = dist->getMax() + 1;

    /* There is no need for shrinking */
    if (delta == 0) {
      return move(dist);
    }

    if (size <= delta) {
      delta = size - 1;
    }

    /* Variables for the distributions ranges */
    int64_t cmin = int64_t(dist->getMin()) - delta;
    int64_t cmax = int64_t(dist->getMax()) - delta;

    /* Determine whether an accumulation is necesary for the shrinking */
    if (cmin < 0) {

      /* Proper shrinked distribution */
      unique_ptr<pmf> clone(new pmf(0, cmax, 1e-4));
      double acum = 0.0;

      /* Perform the accumulation */
      for (int64_t i = cmin; i <= 0; i++) {

        acum += dist->get(i + delta - dist->getMin());

      }

      /* Set the accumulated probability to the first elements of the distribution */
      clone->set(0, acum);

      /* Perform the shifting of the distribution */
      for (int64_t i = 1; i <= cmax; i++) {

        clone->set(i, dist->get(i + delta - dist->getMin()));

      }

      return clone;

    }
    else {

      /* Proper shrinked distribution */
      unique_ptr<pmf> clone(new pmf(cmin, cmax, 1e-4));

      /* Perform the shifting of the distribution */
      for (int64_t i = cmin; i <= cmax; i++) {

        clone->set(i - cmin, dist->get(i + delta - dist->getMin()));

      }

      return clone;

    }

  }

  unique_ptr<pmf> FixedPriorityProbabilitySolver::convolveShrink(unique_ptr<pmf> &first, const unique_ptr<pmf> &second, uint64_t delta_second) {

    /* Size of the distributions starting from 0 */
    uint64_t size1 = first->getMax() + 1;
    uint64_t size2 = second->getMax() + 1;

    if (delta_second > (size1 + size2)) {
      cout << "ERROR" << endl;
      return nullptr;
    }

    /* Shrink the distribution */
    unique_ptr<pmf> first_shrinked = move(convolutionsShrink(first, delta_second));

    /* Convolve the distributions */
    unique_ptr<pmf> clone = move(convolutionsConvolve(first_shrinked, second));

    return clone;

  }

  //void FixedPriorityProbabilitySolver::computeBacklogUntil(unique_ptr<pmf> &starting_distribution, uint64_t id_task, uint64_t until_time) {
  unique_ptr<pmf> FixedPriorityProbabilitySolver::computeBacklogUntil(const unique_ptr<pmf> &starting_distribution, uint64_t id_task, uint64_t until_time) {

    uint64_t previous_time = 0;

    unique_ptr<pmf> start(new pmf(starting_distribution->getMin(), starting_distribution->getMax(), 1e-4));
    for (size_t i = 0; i < starting_distribution->getSize(); i++) {
      start->set(i, starting_distribution->get(i));
    }

    /* Obtain the schedule of activations */
    vector<pair<uint64_t, uint64_t>> sched = task_schedule->getSchedule();

    /* While there is a task and the time is not expired, convolve shrink the next task */
    for (size_t i = 0; (sched[i].second <= until_time) && i < sched.size(); i++) {

      /* The current task has priority or the two tasks have the same priority */
      /* and the release time of the task is smaller than the limit */
      if (task_schedule->hasPriority(sched[i].first, id_task) || 
          (!task_schedule->hasPriority(sched[i].first, id_task) && 
          !task_schedule->hasPriority(id_task, sched[i].first) && 
          (sched[i].second < until_time))) {

        start = move(convolveShrink(start, 
            ((task_schedule->getTaskSet())[sched[i].first])->getComputationTime(), 
            sched[i].second - previous_time));

        previous_time = sched[i].second;

      }

    }

    /* Shrink the distribution of the remaining data */
    start = move(convolutionsShrink(start, until_time - previous_time));

    return start;

  }

  /// @brief Check the conditions to start.
  bool FixedPriorityProbabilitySolver::checkList() {

    /* Check whether there is a linked task */
    if (!linked_flag) {
      EXC_PRINT("ERROR: The solver was called but there is no registered task.");
    }

    /* Check whether the solver was already applied */
    if (solved) {

      /* Present online information */
      if (verbose_flag) {
        cerr << "Solution requested for a problem that has been already solved."
             << endl;
      }

      return false;

    }

    return true;

  }

  /// @brief Prepare the computation.
  void FixedPriorityProbabilitySolver::preProcess() {

    uint64_t r, mr, block_size, M;

    /* Check whether the pre-process was already performed */
    if (pre_process_done) {

      if (verbose_flag) {
        cerr << "The pre-process phase was already performed for the schedule " 
             << task_schedule->getName() << endl;
      }

      return;

    }

    /* Obtain the maximum idle time */
    r = uint64_t(task_schedule->getMaxIdleTime());

    /* Obtain the number of tasks */
    size_t num_tasks = (task_schedule->getTaskSet()).size();

    /* Variable for the PMF of the backlog */
    unique_ptr<pmf> initial_backlog(new pmf(r, r, 1e-4));

    /* Set the PMF of the backlog */
    initial_backlog->set(0, 1.0);

    /* Obtain the (r + 1)-th column of the backlog matrix */
    unique_ptr<pmf> final_backlog = move(computeBacklogUntil(initial_backlog, num_tasks - 1, task_schedule->getHyperperiod()));

    /* Obtain the number of elements in the backlog */
    mr = final_backlog->getSize();

    /* Obtain the size of the submatrices */
    if (r > mr) {
	    block_size = r + 1;
    }
    else {
	    block_size = mr;
    }

    /* Obtain the size of the transition matrix */
    M = 2 * block_size;

    /* Initialize the probability transition matrix */
    MatrixXd P = MatrixXd::Zero(M, M);

    /* Complete the probability transition matrix */
    if (r > mr) {

      /* Construct the shifting columns without cutting the vector */
      for (uint64_t i = 0; i <= (M - block_size); i++) {

        /* Add the repeated shifting block to the matrix */
        P.block(i, r + i, mr, 1) = final_backlog->getElements();

      }

    }
    else {

      /* Construct the shifting columns without cutting the vector */
      for (uint64_t i = 0; i <= (M - block_size); i++) {

        /* Add the repeated shifting block to the matrix */
        P.block(i, r + i, mr, 1) = final_backlog->getElements();

      }

      /* Number of elements to be removed from the pmf resampled vector */
      uint64_t idx = 1;

      /* Construct the columns that shift cutting the pmf resampled vector */
      for (uint64_t i = (M - block_size + 1); i < (M - r); i++) {

        /* Add the cut shifting block to the matrix */
        P.block(i, i + r, (block_size - idx), 1) = (final_backlog->getElements()).block(0, 0, (block_size - idx), 1);

        /* Update the number of rows to cut */
        idx++;

      }

    }

    for (uint64_t i = 0; i < r; i++) {

      /* Variable for the PMF of the backlog */
      unique_ptr<pmf> initial_b_i(new pmf(i, i, 1e-4));

      /* Set the PMF of the backlog */
      initial_b_i->set(0, 1.0);

      /* Obtain the (r + 1)-th column of the backlog matrix */
      unique_ptr<pmf> final_b_i = move(computeBacklogUntil(initial_b_i, num_tasks - 1, task_schedule->getHyperperiod()));

      /* Add the repeated shifting block to the matrix */
      P.block(0, i, final_b_i->getSize(), 1) = final_b_i->getElements();

    }

    /* Initialize submatrices */
    A0 = MatrixXd(block_size, block_size);
    B0 = MatrixXd(block_size, block_size);
    A1 = MatrixXd(block_size, block_size);
    A2 = MatrixXd(block_size, block_size);

    /* Extract submatrices */
    extractSubmatrices(P, block_size, B0, A0, A1, A2);

    /* Initialize matrix */
    R = MatrixXd(A0.rows(), A0.cols());

  }

  /// @brief Post processing after the solution.
  void FixedPriorityProbabilitySolver::postProcess() {

    /* Compute the stationary vector */
    if (!computePi0()) {

      /* Present online information */
      if (task_schedule->getVerboseFlag()) {
        cerr << "WARNING: There were anomalies in the computation of pi0" 
             << endl;
      }

    }

  }

  /// @brief Computes the stationary vector pi0.
  bool FixedPriorityProbabilitySolver::computePi0() {

    /* Check whether there is a linked task */
    if (!task_schedule) {
      EXC_PRINT("ERROR: The computation of pi0 requires a task schedule.");
    }

    /* Check whether the solver has been applied */
    if (!solved) {
      EXC_PRINT_2("ERROR: The solver has not been applied for schedule ", 
          task_schedule->getName());
    }

    /* Check the size of the matrix R */
    if (R.rows() != R.cols()) {
      EXC_PRINT("ERROR: The matrix R has to be square");
    }

    /* The stationary vector was already computed */
    if (post_process_done && verbose_flag) {

      cerr << "WARNING: Computation of pi0 required twice. Schedule " 
           << task_schedule->getName() << endl;

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

    /* MATLAB algorithm */
    MatrixXd temp = (Id - R).inverse();
    RowVectorXd pi = stat(B0 + R * A0);
    pi0 = pi / (pi * temp * u);

    return true;

  }

  /// @brief Fills in the probability map.
  void FixedPriorityProbabilitySolver::fillInProbabilityMap() {

    /* Obtain the number of tasks */
    size_t num_tasks = (task_schedule->getTaskSet()).size();
    
    /* Obtain the schedule of activations */
    vector<pair<uint64_t, uint64_t>> sched = task_schedule->getSchedule();

    uint64_t elems = pi0.cols();

    /* Define the stationary distribution */
    unique_ptr<pmf> stationary(new pmf(0, elems - 1, 1e-4));

    /* Iterate in the range of the stationary distribution */
    for (uint64_t i = 0; i < elems; i++) {

      /* Set the values of the distribution */
      stationary->set(i, pi0(i));

    }

    /* Iterate over the task */
    for (size_t i = 0; i < num_tasks; i++) {

      double prob = 0.0;

      uint64_t el = (((task_schedule->getTaskSet())[i])->getComputationTime())->getMin();
      unique_ptr<pmf> task_resp(new pmf(el, el, 1e-4));

      /* Iterate over the job activations */
      for (size_t j = 0; j < sched.size(); j++) {

        if (sched[j].first == i) {

          /* Compute the backlog until this instance of the task */
          unique_ptr<pmf> job_backlog = move(computeBacklogUntil(stationary, i, sched[j].second));

          unique_ptr<pmf> job_resp = move(getJobResponseTime(job_backlog, i, j));

          task_resp = addDistribution(task_resp, job_resp);

        }

      }

      /* Normalize the sum */
      normalizeDistribution(task_resp);

      /* Shift the distribution */
      shiftDistribution(task_resp, ((task_schedule->getTaskSet())[i])->getActivationTime());

      /* Obtain the deadline probability map */
      DeadlineProbabilityMap *pm = ((task_schedule->getTaskSet())[i])->getProbabilisticDeadlines();
      DeadlineProbabilityMapIter pmi;

      /* Obtain the deadline */
      uint64_t deadline = ((task_schedule->getTaskSet())[i])->getDeadlineStep();

      /* Iterate over the distribution */
      for (uint64_t k = task_resp->getMin(); k <= task_resp->getMax(); k++) {

        /* Accumulate the probability */
        prob += task_resp->get(k - task_resp->getMin());

        /* Set the probability to the multiple of the deadline */
        if ((k % deadline) == 0) {

          /* Obtain the probabilisic deadline */
          if ((pmi = pm->find(k)) != pm->end()) {

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

      for(pmi = pm->begin(); pmi != pm->end(); pmi++) {
        if ((*pmi).second == 0.0) {
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

  }

  /// @brief Executes the iterations of the Latouche method.
  void FixedPriorityProbabilitySolver::applyAlgorithm() {

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

  /// @brief Associates the solver with a task schedule.
  void FixedPriorityProbabilitySolver::registerSchedule(FixedPriorityTaskSchedule *ts) {

    /* Register the task with the probability solver */
    if (!(task_schedule = ts)) {
      EXC_PRINT_2("ERROR: Fixed priority solver used for improper schedule ",
          ts->getName());
    }

    /* Reset the solver */
    reset();

    /* A verbose flag equals to true is not overriden by a false one */
    if (!verbose_flag) {
      verbose_flag = task_schedule->getVerboseFlag();
    }

    /* The solver is linked to a task */
    linked_flag = true;

    /* Present online information */
    if (verbose_flag) {
      cerr << "Schedule " << task_schedule->getName() << " linked to the solver." 
           << endl;
    }

  }

  /// @brief Extract the QBDP submatrices.
  void FixedPriorityProbabilitySolver::extractSubmatrices(const MatrixXd &P,
      uint64_t dim, MatrixXd &B, MatrixXd &A0, MatrixXd &A1, MatrixXd &A2) {

    /* The obtained matrix must be transposed due to the definition in the original paper */
    MatrixXd T = P.transpose();

    B = T.block(0, 0, dim, dim);
    A0 = T.block(dim, 0, dim, dim);
    A1 = T.block(dim, dim, dim, dim);
    A2 = T.block(0, dim, dim, dim);

  }

}
