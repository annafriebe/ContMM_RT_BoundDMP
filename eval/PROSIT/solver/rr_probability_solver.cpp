/*!
 * @file    rr_probability_solver.cpp
 * 
 * @brief   This class contains the implementation for generic probability
 *          solvers for resource reservation tasks. In essence, a probability
 *          solver wraps the different algorithms used to compute probabilistic
 *          deadlines. The solver for resource reservation is one to one linked
 *          with a task.
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
#include "rr_probability_solver.hpp"

using namespace PrositAux;
using namespace std;

namespace PrositCore {

  /// @brief Check the conditions to start.
  bool ResourceReservationProbabilitySolver::checkList() {

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

    /* Check whether the deadline step has been defined */
    if (task_descriptor->getDeadlineStep() == 0) {
      EXC_PRINT_2("ERROR: The solver has been called, however the deadline "
          "step has not been defined for task ", task_descriptor->getName());
    }

    /* Check whether the probabilistic deadlines has been inserted */
    if (task_descriptor->getProbabilisticDeadlines()->empty()) {
      EXC_PRINT_2("ERROR: The solver has been called, however the "
          "probabilistic deadlines has not been set for task ", 
          task_descriptor->getName());
    }

    return true;

  }

  /// @brief Associates the solver with a task descriptor.
  void ResourceReservationProbabilitySolver::registerTask(ResourceReservationTaskDescriptor *td) {

    /* Register the task with the probability solver */
    if (!(task_descriptor = td)) {
      EXC_PRINT_2("ERROR: Resource reservation solver used for improper task ",
          td->getName());
    }

    /* Reset the solver */
    reset();

    /* Properly set the verbose flag */
    verbose_flag = task_descriptor->getVerboseFlag();

    /* The solver is linked to a task */
    linked_flag = true;

    /* Present online information */
    if (verbose_flag) {
      cerr << "Task " << task_descriptor->getName() << " linked to the solver." 
           << endl;
    }

  }

}
