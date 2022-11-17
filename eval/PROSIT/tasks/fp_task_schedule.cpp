/*!
 * @file    fp_task_schedule.cpp
 * 
 * @brief   
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
#include "../solver/fp_probability_solver.hpp"
#include "fp_task_schedule.hpp"

using namespace PrositAux;
using namespace std;

namespace PrositCore {

  /// @brief Constructor.
  FixedPriorityTaskSchedule::FixedPriorityTaskSchedule(string nm, FixedPriorityTaskSet tasksetd) {

    /* Define the name */
    name = nm;

    /* Define the verbose flag */
    verbose_flag = false;

    /* Define the taskset */
    taskset = move(tasksetd);

    /* Compute the hyperperiod of the task set */
    computeHyperperiod();

    /* Obtain the schedule of jobs */
    buildSchedule();

    /* Compute the maximum idle time */
    computeMaximumIdleTime();

  }

  /*FixedPriorityTaskSchedule::FixedPriorityTaskSchedule(const FixedPriorityTaskSchedule &s) {
          name = s.name;
          verbose_flag = s.verbose_flag;
          solved = s.solved;
          hyperperiod = s.hyperperiod;
          maximum_idle_time = s.maximum_idle_time;
          schedule = s.schedule;
          taskset = s.taskset;
          probability_solver = new FixedPriorityProbabilitySolver(*s.probability_solver);
  }*/

  /// @brief Destructor.
  FixedPriorityTaskSchedule::~FixedPriorityTaskSchedule() = default;

  /// @brief Obtain the probability associated with a deadline.
  double FixedPriorityTaskSchedule::getProbability(DeadlineUnit deadline, uint64_t task_id) {

    /* Check if the task is associated to a solver */
    if (!probability_solver) {
      EXC_PRINT_2("ERROR: No solver linked to task", name);
    }

    /* Obtain the iterator to the deadline */
    DeadlineProbabilityMapIter it = (taskset[task_id]->getProbabilisticDeadlines())->find(deadline);

    /* Check if the task has been solved */
    if (!solved) {

      if (verbose_flag) {
        cout << "Probability requested for unsolved task. Now solving..." << endl;
      }

      /* Compute the probability of respecting the deadline */
      computeProbability();

    }

    /* Check if the given deadline belong to the map */
    if (it == (taskset[task_id]->getProbabilisticDeadlines())->end()) {
      EXC_PRINT_2("ERROR: Deadline does not exist for task ", name);
    }

    /* Return the probability */
    return it->second;

  }

  /// @brief Set the probability solver.
  void FixedPriorityTaskSchedule::setSolver(unique_ptr<FixedPriorityProbabilitySolver> psd) {

    /* Set the probability solver */
    probability_solver = move(psd);

    /* Declare the schedule as unsolved */
    solved = false;

    /* Register the task in the solver */
    probability_solver->registerSchedule(this);

  }

  /// @brief Compute the least common multiple of two input numbers.
  uint64_t FixedPriorityTaskSchedule::computeLeastCommonMultiple(uint64_t n1, uint64_t n2) {

    uint64_t lcm;

    /* Determine the greater input argument */
    lcm = (n1 > n2) ? n1 : n2;

    /* The computed LCM is not multiple of both input arguments */
    while (((lcm % n1) != 0) || ((lcm % n2) != 0)) {

      /* Increase the least common multiple */
      lcm++;

    }

    return lcm;

  }

  /// @brief Compute the hyperperiod of a taskset.
  void FixedPriorityTaskSchedule::computeHyperperiod() {

    /* Obtain the number of tasks */
    size_t num_tasks = taskset.size();

    /* There is only 1 task in the set */
    if (num_tasks == 1) {

      /* The hyperperiod is equal to the period of the task */
      hyperperiod = taskset[0]->getPeriod();

    }
    else {

      /* Iterate over the tasks in the taskset */
      for (size_t i = 1; i < num_tasks; i++) {

        if (i == 1) {

          /* Compute the hyperperiod between the first 2 tasks */
          hyperperiod = computeLeastCommonMultiple(taskset[0]->getPeriod(), taskset[1]->getPeriod());

        }
        else {

          /* Compute the hyperperiod between the current hiperperiod and the next task */
          hyperperiod = computeLeastCommonMultiple(hyperperiod, taskset[i]->getPeriod());

        }

      }

    }

  }

  /// @brief Obtain the stream of activations.
  void FixedPriorityTaskSchedule::buildSchedule() {

    /* Obtain the number of tasks */
    size_t num_tasks = taskset.size();

    /* Define the initial time instant and task */
    uint64_t current_time = 0;
    uint64_t current_task = 0;
    bool proceed = true;

    /* There are jobs activations to include */
    while (proceed) {

      bool go_ahead = true;
      uint64_t task_copy = current_task;

      /* Iterate over the time */
      for (uint64_t i = current_time; go_ahead && i < (hyperperiod - 1); i++) {

        /* Iterate over the tasks */
        for (uint64_t j = task_copy; go_ahead && j < num_tasks; j++) {

          /* A job activation has occurred */
          if (((i % taskset[j]->getPeriod()) - taskset[j]->getActivationTime()) == 0) {

		        /* Store the task and time of the job activation */
            go_ahead = false;
		        current_task = j;
            current_time = i;

          }

        }

        task_copy = 0;

      }

      /* The hyperperiod has expired */
      if (go_ahead) {

        proceed = false;

      }
      else {

        /* Insert the job in the schedule */
        schedule.push_back(make_pair(current_task, current_time));

      }

	    /* Increase the time and restart the task */
      if (num_tasks == (current_task + 1)) {

	      current_time++;
	      current_task = 0;

	    }
      else {

	      current_task++;

	    }

    }

  }

  /// @brief Determine whether a task has priority over another task.
  bool FixedPriorityTaskSchedule::hasPriority(uint64_t task1, uint64_t over_task2) {

    if (taskset[task1]->getDeadlineStep() < taskset[over_task2]->getDeadlineStep()) {
	    return true;
    }
    else {
	    return false;
    }

  }

  /// @brief Compute the maximum idle time in any hiperperiod.
  void FixedPriorityTaskSchedule::computeMaximumIdleTime() {

    /* Initialize the variables */
    int64_t prev_time = schedule[0].second;
    int64_t prev_duration = (taskset[schedule[0].first]->getComputationTime())->getMin();
    int64_t W_min = prev_duration;
    int64_t C_min = 0;

    /* Iterate over the schedule */
    for (size_t i = 0; i < schedule.size(); i++) {

      /* The accumulated backlog is updated */
      W_min = (W_min > int64_t(schedule[i].second) - prev_time) ? W_min - (int64_t(schedule[i].second) - prev_time) : 0;

      /* Update the release time and computation time of the job activation */
      prev_duration = (taskset[schedule[i].first]->getComputationTime())->getMin();
      prev_time = schedule[i].second;

      /* Accumulate the backlog of the activated jobs */
      W_min += prev_duration;

      /* Accumulate the minimum computation time of the activated jobs */
      C_min += (taskset[schedule[i].first]->getComputationTime())->getMin();
      
    }

    /* The accumulated backlog is updated */
    W_min = (W_min > int64_t(hyperperiod) - prev_time) ? W_min - (int64_t(hyperperiod) - prev_time) : 0;

    /* Compute the maximum idle time in the hyperperiod */
    maximum_idle_time = int64_t(hyperperiod) + W_min - C_min;

  }

  /// @brief Compute the probability of respecting the deadlines.
  void FixedPriorityTaskSchedule::computeProbability() {

    /* Check if the probability solver has been set */
    if (!probability_solver) {
      EXC_PRINT_2("ERROR: The probability solver has not been set for "
          "task ", name);
    }

    /* Check if the task is already solver */
    if (probability_solver->is_solved()) {

      /* Update the solved flag */
      solved = true;

      return;

    }

    /* Call the probability solver */
    probability_solver->solve();

  }

}
