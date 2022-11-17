/*!
 * @file    rr_task_descriptor.cpp
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
#include "rr_task_descriptor.hpp"

using namespace PrositAux;
using namespace Eigen;
using namespace std;

namespace PrositCore {

  /// @brief Constructor.
  ResourceReservationTaskDescriptor::ResourceReservationTaskDescriptor(std::string nm, 
      int num_modesd,
      std::vector<std::unique_ptr<PrositAux::pmf>> Cd, 
      std::unique_ptr<PrositAux::pmf> Zd, 
      Eigen::MatrixXd transition_matrixd,
      const uint64_t Qsd, 
      const uint64_t Tsd, 
      uint64_t granularityd,
      uint64_t deadlined,
      std::string algorithm) : 
      GenericTaskDescriptor(nm, std::move(Zd), deadlined, algorithm), 
      solved(false), 
      num_modes(num_modesd),
      C(move(Cd)), 
      transition_matrix(transition_matrixd),
      Qs(Qsd),
      Ts(Tsd),
      granularity(granularityd),
      Qmin(0),
      Qmax(Tsd) {

    /* Check the correct value for the bandwidth */
    if ((double(Qs) / double(Ts)) > 1.0) {
      EXC_PRINT_2("ERROR: The server period is too small for task ", name);
    }

    /* Update the bandwidth */
    bandwidth = double(Qs) / double(Ts);

    /* Check whether the task is periodic */
    if (Z->getMin() == Z->getMax()) {
      is_periodic = true;
      period = Z->getMin() * Ts;
    }

    /* Perform the padding to the distributions */
    if (num_modes > 1) {

      /* Variables for the distributions ranges */
      uint64_t cmin = numeric_limits<uint64_t>::max();
      uint64_t cmax = numeric_limits<uint64_t>::min();

      /* Iterate over the number of modes */
      for (int i = 0; i < num_modes; i++) {

        /* Obtain the global minimum */
        if (C[i]->getMin() < cmin) {
          cmin = C[i]->getMin();
        }

        /* Obtain the global maximum */
        if (C[i]->getMax() > cmax) {
          cmax = C[i]->getMax();
        }

      }

      /* Iterate over the number of modes */
      for (int i = 0; i < num_modes; i++) {

        /* Variable for the new vector of elements */
        VectorXd tmp_elem = VectorXd::Zero(cmax - cmin + 1);

        /* Build the new vector of elements */
        tmp_elem.block(C[i]->getMin() - cmin, 0, C[i]->getSize(), 1) = C[i]->getElements();

        /* Set the new elements of the distribution */
        C[i]->setMin(cmin);
        C[i]->setMax(cmax);
        C[i]->setSize(cmax - cmin + 1);
        C[i]->setElements(tmp_elem);

      }

    }

    /* Initialize the bounds */
    bounds_inited = false;

  }

  /// @brief Destructor.
  ResourceReservationTaskDescriptor::~ResourceReservationTaskDescriptor() = default;

  /// @brief Set the probability solver.
  void ResourceReservationTaskDescriptor::setSolver(unique_ptr<ResourceReservationProbabilitySolver> psd) {

    /* Set the probability solver */
    probability_solver = move(psd);

    /* Declare the task as unsolved */
    solved = false;

    /* Register the task in the solver */
    probability_solver->registerTask(this);

  }

  void ResourceReservationTaskDescriptor::resetSolver() {

    solved = false;

    /* Reset the solver to ensure the new calculation of the probability */
    probability_solver->reset();

    /* Set to zero the probabilities */
    cleanProbabilisticDeadlines();

  }

  /// @brief Compute the probability of respecting the deadlines.
  void ResourceReservationTaskDescriptor::computeProbability() {

    /* Check if the probability solver has been set */
    if (!probability_solver) {
      EXC_PRINT_2("ERROR: The probability solver has not been set for "
          "task ", name);
    }

    /* Check if the probabilistic deadline map is empty */
    if (probabilistic_deadlines.empty()) {
      EXC_PRINT_2("ERROR: No deadline specified for task ", name);
    }

    /* Check if the task is already solver */
    if (probability_solver->is_solved()) {

      /* Update the solved flag */
      solved = true;

      return;

    }

    /* Call the probability solver */
    probability_solver->solve();

//        DeadlineProbabilityMapIter it = probabilistic_deadlines.begin();

        /* Iterate over the map */
//        while (it != probabilistic_deadlines.end()) {

          /* Set the probability to zero */
//          cout << it->first << " - " << it->second << endl;

          /* Increment the ierator to point to next entry */
//          it++;

//        }

  }

  /// @brief Set the deadline step.
  void ResourceReservationTaskDescriptor::setDeadlineStep(DeadlineUnit ds) {

    /* Check if the deadline step is a multiple of the server period */
    if ((ds != 0) && ((ds % Ts) > 1)) {
      EXC_PRINT_2("ERROR: The deadline step has to be a multiple of the server "
          "period for ", name);

    }

    /* Set the reservation period */
    GenericTaskDescriptor::setDeadlineStep(ds);

  }

  /// @brief Obtain the probability associated with a deadline.
  double ResourceReservationTaskDescriptor::getProbability(DeadlineUnit deadline) {

    /* Check if the task is associated to a solver */
    if (!probability_solver) {
      EXC_PRINT_2("ERROR: No solver linked to task", name);
    }

    /* Obtain the iterator to the deadline */
    DeadlineProbabilityMapIter it = probabilistic_deadlines.find(deadline);

    /* Check if the task has been solved */
    if (!solved) {

      if (verbose_flag) {
        cout << "Probability requested for unsolved task. Now solving..." << endl;
      }

      /* Compute the probability of respecting the deadline */
      computeProbability();

    }

    /* Check if the given deadline belong to the map */
    if (it == probabilistic_deadlines.end()) {
      EXC_PRINT_2("ERROR: Deadline does not exist for task ", name);
    }

    /* Return the probability */
    return it->second;

  }

  /// @brief Evaluate the QoS function for a given budget.
  double ResourceReservationTaskDescriptor::QoSEvaluation(DeadlineUnit Q) {

    /* Reset the solver to ensure the new calculation of the probability */
    probability_solver->reset();

    double prob = 0.0;
    
    /* Compute the probability of respecting the deadline */
    if (Q != 0) {

      /* Save the previous value of the budget */
      uint64_t prev_Q = Qs;
      uint64_t prev_granularity = 1;

      /* Update the value of the budget */
      setBudget(Q);

      /* Update the granularity to be the budget */
      if (strcmp(algorithm.c_str(), "analytic") == 0) {
        prev_granularity = granularity;
        granularity = Qs;
      }

      /* Compute the probability with the new budget */
      prob = getProbability(period);

      /* Restore the previous budget */
      setBudget(prev_Q);

      /* Restore the granularity */
      if (strcmp(algorithm.c_str(), "analytic") == 0) {
        granularity = prev_granularity;
      }

    }

    /* Evaluate the QoS function with the new probability */
    return q->evaluateQuality(prob);

  }

  bool ResourceReservationTaskDescriptor::inverseQoSEvaluation(double p, uint64_t &Q, bool ceil) {

    if (QoSEvaluation(Ts) < p) {

      if (ceil) {
        cerr << "WARNING: Unfeasible to fullfill the QoS target (" << p << ") for task " << name << endl;
        return false;
      }
      else {
        Q = Ts;
        return true;
      }

    }

    if (QoSEvaluation(0) > p) {

      if (!ceil) {
        cerr << "WARNING: Unfeasible to fullfill the QoS target (" << p << ") for task " << name << endl;
        return false;
      }
      else {
        Q = 0;
        //cout << "ECCOMI" << endl;
        return false;
      }

    }

    uint64_t Qm = 0, QM = Ts;

    while (QM > (Qm + 1)) {

      //cout << "Qm: " << Qm << " - QM: " << QM << endl;
      if (QoSEvaluation((Qm + QM) / 2) < (p - 1e-8))
        Qm = (Qm + QM) / 2;
      else
        QM = (QM + Qm) / 2;

    }

    if (ceil) {
      Q = QM;
    }
    else {
      if (QoSEvaluation(QM) == p)
        Q = QM;
      else
        Q = Qm;
    }
    //cout << "Budget: " << Q << endl;

    return true;

  }

  bool ResourceReservationTaskDescriptor::identifyBounds(double q_min, double q_max) {

    bool res1, res2;

    if (q_min > q_max)
      EXC_PRINT_2("Bad bounds made identify_bounds fail for task", name);

    res1 = inverseQoSEvaluation(q_min, Qmin, true);
    res2 = inverseQoSEvaluation(q_max, Qmax, false);

    if (Qmin > Qmax) {
      return false;
    }

    QoSmin = QoSEvaluation(Qmin);

    QoSmax = QoSEvaluation(Qmax);

    bounds_inited = true;

    return res1 && res2;

  }

  void ResourceReservationTaskDescriptor::display(string, GenericTaskDescriptor* td, const vector<double> &probability, const vector<double> &quality, const vector<long long> &time, int index) {

    ResourceReservationTaskDescriptor *t;

    if (!(t = dynamic_cast<ResourceReservationTaskDescriptor *>(td))) {
      EXC_PRINT_2("Impossible to cast task GenericTaskDescriptor to ResourceReservationTaskDescriptor for task ", td->getName());
    }

    /* Present the results */
    cout << setw(8) << t->getName() 
         << setw(18) << t->getBudget() << fixed
         << setw(19) << setprecision(2) << double(t->getBudget()) / double(t->getServerPeriod()) 
         << setw(35) << setprecision(20) << probability[index]
         << setw(15) << setprecision(4) << quality[index] 
         << setw(15) << time[index] << " us" << endl;

    t->inf_norm = min<double>(quality[index], t->inf_norm);

  }

}
