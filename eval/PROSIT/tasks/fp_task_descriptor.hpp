/*!
 * @file    fp_task_descriptor.hpp
 * 
 * @brief   This class defines the header for the fixed priority task 
 *          descriptor, where the tasks are scheduled using fixed priorities.
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
#ifndef FP_TASK_DESCRIPTOR_HPP
#define FP_TASK_DESCRIPTOR_HPP

#include "generic_task_descriptor.hpp"

namespace PrositCore {

  class FixedPriorityTaskDescriptor : public GenericTaskDescriptor {

    protected:

      ///< Distribution of the computation time.
      std::unique_ptr<PrositAux::pmf> C;

      uint64_t priority;           ///< Scheduling priority of the task.
      uint64_t activation_time;    ///< Activation time of the task.

    public:
  
      /// @brief Aperiodic constructor.
      ///
      /// This is the constructor for aperiodic fixed priority tasks.
      ///
      /// @param nm is the unique identifier for the task.
      /// @param Cd is the pointer to the distribution of the computation time.
      /// @param Zd is the pointer to the distribution of the interarrival time.
      /// @param priorityd is the scheduling priority (in the range 1..99).
      /// @param algorithm is the solver to be set.
      FixedPriorityTaskDescriptor(std::string nm, 
          std::unique_ptr<PrositAux::pmf> Cd, 
          std::unique_ptr<PrositAux::pmf> Zd, 
          uint64_t priorityd, 
          uint64_t deadlined, 
          uint64_t activationd) : 
          GenericTaskDescriptor(nm, std::move(Zd), deadlined, "logarithmic"),
          C(move(Cd)),
          activation_time(activationd) {

        /* Check the correct range of the priority */
        if ((priorityd < 1) || (priorityd > 99)) {
          EXC_PRINT_2("ERROR: Priority out of range for task ", nm);
        }

        /* Set the priority */
        priority = priorityd;

        /* Set the period */
        period = Z->getMin();
        is_periodic = true;

        /* Set the deadline */
        deadline_step = deadlined;

      }

      /// @brief Destructor.
      ///
      /// This is the destructor of the fixed priority task descriptor.
      virtual ~FixedPriorityTaskDescriptor() { }

      /// @brief Obtain a copy of the computation time distribution.
      ///
      /// This function is used to return a copy of the PMF of the computation 
      /// time.
      ///
      /// @return a copy of the pmf related to the computation time.
      std::unique_ptr<PrositAux::pmf> const& getComputationTime() const {
        return C;
      }

      /// @brief Set the task priority.
      ///
      /// This function sets a new priority for the task.
      ///
      /// @param priority is the new task priority.
      ///
      /// @return old_prio is the old task priority.
      uint64_t setPriority(uint64_t priorityd) {

        /* Store the old task priority */
        unsigned int old_prio = priority;

        /* Set the new priority */
        priority = priorityd;

        return old_prio;

      }

      /// @brief Obtain the task priority.
      ///
      /// This function returns the priority of the task.
      ///
      /// @return priority is the task priority.
      uint64_t getPriority() const {
        return priority;
      }

      /// @brief Obtain the activation time of the task.
      ///
      /// This function returns the activation time of the task.
      ///
      /// @return activation_time is the activation time of the task.
      uint64_t getActivationTime() const {
        return activation_time;
      }

      int getNumModes() const {
        return 1;
      }

      /// @brief Displays results for a single task.
      ///
      /// @param td is the ...
      /// @param probability is the ...
      /// @param quality is the ...
      /// @param time is the ...
      /// @param show_time is the ...
      /// @param index is the ...
      virtual void display(std::string name_taskset,
          GenericTaskDescriptor* td, 
          const std::vector<double> &probability, 
          const std::vector<double> &quality, 
          const std::vector<long long> &time, 
          int index);

  };

}

#endif
