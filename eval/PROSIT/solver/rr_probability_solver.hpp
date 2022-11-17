/*!
 * @file    rr_probability_solver.hpp
 * 
 * @brief   This class contains the abstract definitions for generic 
 *          probability solvers for resource reservation tasks. In 
 *          essence, a probability solver wraps the different algorithms 
 *          used to compute probabilistic deadlines. The solver for resource 
 *          reservation is one to one linked with a task.
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
#ifndef RR_PROBABILITY_SOLVER_HPP
#define RR_PROBABILITY_SOLVER_HPP

#include "../tasks/rr_task_descriptor.hpp"
#include "probability_solver.hpp"

namespace PrositCore {

  /* Forward declaration */
  class ResourceReservationTaskDescriptor;

  class ResourceReservationProbabilitySolver : public ProbabilitySolver {

    protected:

      ///< Descriptor of only task managed by this type of solver.
      ResourceReservationTaskDescriptor *task_descriptor; 

      /// @brief Check the conditions to start.
      ///
      /// This function determines whether everything is ok to start the 
      /// pre-process phase.
      bool checkList();

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for resource reservation probability solvers.
      ResourceReservationProbabilitySolver() : 
          ProbabilitySolver(), 
          task_descriptor(0) { }

      /// @brief Destructor.
      ///
      /// This is the destructor for resource reservation probability solvers.
      virtual ~ResourceReservationProbabilitySolver() { }

      /// @brief Register the task to the solver.
      ///
      /// This function add the task to the set of tasks the solver applies 
      /// to. A resource reservation solver applies to a resource reservation
      /// task descriptor. The call to this method resets previous computations
      /// and sets the verbose flag. An exception is thrown if the task
      /// descriptor does not correspond to a resource reservation task.
      ///
      /// @param td descriptor of the task.
      virtual void registerTask(ResourceReservationTaskDescriptor *td);

  };

}
#endif
