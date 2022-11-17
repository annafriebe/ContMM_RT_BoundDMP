/*!
 * @file    fp_task_schedule.hpp
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
#ifndef FP_TASK_SCHEDULE_HPP
#define FP_TASK_SCHEDULE_HPP

#include "../tasks/fp_task_descriptor.hpp"

namespace PrositCore {

  /* Forward declaration */
  class FixedPriorityProbabilitySolver;

  class FixedPriorityTaskSchedule {

    protected:

      std::string name;             ///< Name of the schedule.
      bool verbose_flag;            ///< Flag to print out online information.
      bool solved;                  ///< The solver has been called.

      uint64_t hyperperiod;         ///< Hyperperiod of the taskset.
      int64_t maximum_idle_time;    ///< Maximum possible value of the idle time in any hyperperiod.

      ///< Vector of jobs activations.
      std::vector<std::pair<uint64_t, uint64_t>> schedule;

      ///< Vector of fixed priority tasks.
      FixedPriorityTaskSet taskset;

      ///< Pointer to the object containing the solution algorithm for 
      ///< probabilities.
      std::unique_ptr<FixedPriorityProbabilitySolver> probability_solver;

      /// @brief Compute the least common multiple of two input numbers.
      ///
      /// This function computes the least common multiple of 2 given numbers.
      /// The least common multiple is the smallest positive integer that is 
      /// divisible by the two input numbers.
      ///
      /// @param a is the first number.
      /// @param b is the second number.
      ///
      /// @return lcm is the least common multiple of the input arguments.
      uint64_t computeLeastCommonMultiple(uint64_t n1, uint64_t n2);

      /// @brief Compute the hyperperiod of a taskset.
      ///
      /// This function computes the hyperperiod of a taskset. The hyperperiod
      /// is defined as the least common multiple of the periods of all the 
      /// periodic tasks in a taskset.
      void computeHyperperiod();

      /// @brief Obtain the stream of activations.
      ///
      /// This function generates an ordered sequence with the jobs activated 
      /// during a hyperperiod. It produces a vector of pairs containing the 
      /// identifier of the task and the activation time of the job.
      void buildSchedule();

      /// @brief Compute the maximum idle time in any hiperperiod.
      ///
      /// This function computes the maximum possible value of the idle
      /// time in any hyperperiod.
      void computeMaximumIdleTime();

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for fixed-priority task schedule.
      ///
      /// @param tasksetd is the taskset of fixed-priority tasks.
      FixedPriorityTaskSchedule(std::string nm, FixedPriorityTaskSet tasksetd);

      /// @brief Constructor.
      ///
      /// This is the copy constructor of the tasks schedule.
      ///
      /// @param s is an already existing schedule.
      //FixedPriorityTaskSchedule(const FixedPriorityTaskSchedule &s);

      /// @brief Destructor.
      ///
      /// This is the destructor of the resource reservation task descriptor.
      virtual ~FixedPriorityTaskSchedule();

      /// @brief Obtain the name of the task.
      ///
      /// This function returns the name of the task.
      ///
      /// @return a string containing the name of the task.
      std::string getName() const {
        return name;
      }

      /// @brief Set the verbose flag to a specified value.
      ///
      /// This function sets the verbose flag to the value specified as 
      /// parameter.
      ///
      /// @param verbose_flagd is the desidered value for the verbose flag.
      ///
      /// @return the old value of the verbose flag.
      bool setVerboseFlag(bool verbose_flagd) {

        /* Store the old verbose flag */
        bool current = verbose_flag;

        /* Set the new verbose flag */
        verbose_flag = verbose_flagd;

        return current;

      }

      /// @brief Obtain the verbose flag.
      ///
      /// This function returns the current verbose flag.
      ///
      /// @return the current value of the verbose flag.
      bool getVerboseFlag() const {
        return verbose_flag;
      }

      /// @brief Set the probability solver.
      ///
      /// This function sets the external object that computes the probability.
      /// The probability computation has to be made anew. The ownership of 
      /// the solver is not taken by the task descriptor. If the solver is not 
      /// appropriate for the task (e.g., if we use a FP solver for a resource 
      /// reservation task), an exception is thrown by the solver during the 
      /// task registration phase.
      ///
      /// @param psd is a pointer to the probability solver that has to be used.
      void setSolver(std::unique_ptr<FixedPriorityProbabilitySolver> psd);

      /// @brief Obtain the hyperperiod of the task set.
      ///
      /// This function returns the hyperperiod of the taskset.
      ///
      /// @return the hyperperiod of the taskset.
      uint64_t getHyperperiod() const {
        return hyperperiod;
      }

      /// @brief Compute the maximum idle time in any hiperperiod.
      ///
      /// This function returns the maximum amount of time that the processor
      /// is idle while serving a set of tasks.
      ///
      /// @return the maximum idle time of the taskset.
      int64_t getMaxIdleTime() const {
        return maximum_idle_time;
      }

      /// @brief Obtain the schedule of job activations.
      ///
      /// This function returns the sequence of jobs activated within a 
      /// hyperperiod.
      ///
      /// @return the schedule of the taskset.
      std::vector<std::pair<uint64_t, uint64_t>> getSchedule() const {
        return schedule;
      }

      /// @brief Obtain a copy of the computation time distribution.
      ///
      /// This function is used to return a copy of the PMF of the computation 
      /// time.
      ///
      /// @return a copy of the pmf related to the computation time.
      FixedPriorityTaskSet getTaskSet() {
        return taskset;
      }

      /// @brief Determine whether a task has priority over another task.
      ///
      /// This function determines whether a job of a given task has priority
      /// over a job of another given task.
      ///
      /// @param task1 is the tasks of the first job.
      /// @param task2 is the task of the second job.
      ///
      /// @return true if the job of task1 has priority over the job of task2.
      bool hasPriority(uint64_t task1, uint64_t over_task2);

      /// @brief Obtain the probability associated with a deadline.
      ///
      /// This function returns the probability associated with a given 
      /// deadline. The function computeProbability() is implicitly called 
      /// if not called before. An exception is thrown if the deadline has
      /// not been registered.
      ///
      /// @param deadline is the deadline for which the probability is required.
      ///
      /// @return the requested probability.
      double getProbability(DeadlineUnit deadline, uint64_t task_id);

      /// @brief Compute the probability of respecting the deadlines.
      ///
      /// These probabilities are computed for a given configuration of the
      /// scheduling parameters. The deadlines have been previously registered
      /// by calling the function insertDeadline(). An exception is thrown if 
      /// the solver is not set.
      void computeProbability();

  };

}

#endif
