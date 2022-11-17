/*!
 * @file    generic_task_descriptor.cpp
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
#include "generic_task_descriptor.hpp"

using namespace PrositAux;
using namespace std;

namespace PrositCore {

  /// @brief Set the deadline step.
  void GenericTaskDescriptor::setDeadlineStep(DeadlineUnit ds) {

    /* Check if the task is already defined */
    if ((deadline_step != 0) && (!probabilistic_deadlines.empty())) {
      EXC_PRINT_2("ERROR: Deadline step reset in presence of probabilistic "
          "deadlines for task ", name);
    }

    /* Set the deadline step */
    deadline_step = ds;

    return;

  }

  /// @brief Define a deadline in the map.
  void GenericTaskDescriptor::insertDeadline(DeadlineUnit deadline) {

    /* The deadline step has not been defined */
    if (deadline_step == 0) {
      EXC_PRINT_2("ERROR: Deadline inserted before defining the step for "
          "task ", name);
    }

    /* Create the deadline entry */
    pair<DeadlineUnit, double> entry(deadline, 0.0);

    /* Insert the deadline in the map */
    if (!(probabilistic_deadlines.insert(entry).second)) {
      EXC_PRINT_2("ERROR: Cannot create deadline for task ", name);
    }

    return;

  }

}
