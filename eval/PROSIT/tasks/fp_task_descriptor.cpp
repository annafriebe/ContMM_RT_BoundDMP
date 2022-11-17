/*!
 * @file    fp_task_descriptor.cpp
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
#include "fp_task_descriptor.hpp"

using namespace PrositAux;
using namespace std;

namespace PrositCore {

  void FixedPriorityTaskDescriptor::display(string name_taskset, GenericTaskDescriptor* td, const vector<double> &probability, const vector<double> &quality, const vector<long long> &time, int index) {

    FixedPriorityTaskDescriptor *t;

    if (!(t = dynamic_cast<FixedPriorityTaskDescriptor *>(td))) {
      EXC_PRINT_2("Impossible to cast task GenericTaskDescriptor to FixedPriorityTaskDescriptor for task ", td->getName());
    }

    /* Present the results */
    cout << setw(11) << name_taskset 
         << setw(9) << t->getName() 
         << setw(10) << t->getPriority()
         << setw(14) << t->getPeriod()
         << setw(13) << t->getDeadlineStep() << fixed 
         << setw(28) << setprecision(20) << probability[index]
         << setw(12) << setprecision(4) << quality[index] 
         << setw(13) << time[index] << " us" << endl;

    t->inf_norm = min<double>(quality[index], t->inf_norm);

  }

}
