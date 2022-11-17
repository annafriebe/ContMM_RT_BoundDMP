/*!
 * @file    generic_task_descriptor.hpp
 * 
 * @brief   This class contains the root of the hierarchy defining the task
 *          descriptor: the GenericTaskDescriptor, which models both periodic 
 *          and aperiodic tasks depending on a flag.
 *          This descriptor only contains the timing information. The hierarchy 
 *          is specialised on the basis of: scheduling algorithm and solver 
 *          family.
 *          The key method in compute_probability(), which computes the 
 *          probability of respecting a sequence of deadlines that have to be
 *          registered before calling the method. This sequence is composed 
 *          of multiples of a fixed basic deadline.
 *          The computation of the probability relies on an externally provided 
 *          solver (see the documentation of the PrositCore::Solver class).
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
#ifndef GENERIC_TASK_DESCRIPTOR_HPP
#define GENERIC_TASK_DESCRIPTOR_HPP

#include "../qos/quadratic_qos_function.hpp"
#include "../qos/linear_qos_function.hpp"
#include "../utils/prosit_types.hpp"
#include "../utils/pmf.hpp"

namespace PrositCore {

  /* Forward declaration */
  class QoSFunction;

  /* Forward declaration */
  class ResourceReservationProbabilitySolver;

  class GenericTaskDescriptor {

    protected:

      std::string name;   ///< Name of the task.

      ///< Distribution of the interarrival time.
      std::unique_ptr<PrositAux::pmf> Z;

      bool verbose_flag;            ///< Flag to print out online information.
      bool is_periodic;             ///< Flag for periodic tasks.
      uint64_t period;              ///< Period (only for periodic tasks).
      uint64_t deadline;            ///< Deadline.
      DeadlineUnit deadline_step;   ///< Step for the probabilistic deadlines.

      ///< Map associating deadlines with probabilities.
      DeadlineProbabilityMap probabilistic_deadlines;

    public:

      double inf_norm = 1e38;    ///< Used to calculate the QoS.
      std::string algorithm;     ///< Used to set the solver.

      std::unique_ptr<QoSFunction> q;    ///< Pointer to the QoS function.
      std::string qos_type;              ///< Type of QoS function..

      /// @brief Constructor.
      ///
      /// This is the constructor for generic fixed-priority tasks.
      ///
      /// @param nm is the unique identifier for the task.
      /// @param Cd is the vector of pointers to the distribution of the computation time.
      ///        The descriptor takes ownership of this pointer.
      /// @param Zd is the pointer to the distribution of the interarrival time.
      ///        The descriptor takes ownership of this pointer.
      /// @param deadlined is the deadline of the task.
      GenericTaskDescriptor(std::string nm, 
          std::unique_ptr<PrositAux::pmf> Zd,
          uint64_t deadlined, 
          std::string algorithm): 
          name(nm), 
          Z(move(Zd)), 
          verbose_flag(false), 
          is_periodic(false), 
          period(0), 
          deadline(deadlined), 
          deadline_step(0), 
          algorithm(algorithm) { }

      /// @brief Destructor.
      ///
      /// This is the destructor of the task descriptor, which is virtual being
      /// the class polymorphic.
      virtual ~GenericTaskDescriptor() { }

      /// @brief Obtain the name of the task.
      ///
      /// This function returns the name of the task.
      ///
      /// @return a string containing the name of the task.
      std::string getName() const {
        return name;
      }

      virtual int getNumModes() const = 0;
      virtual void setSolver(std::unique_ptr<ResourceReservationProbabilitySolver>) { }

      /// @brief Obtain a copy of the interarrival time distribution.
      ///
      /// This function is used to return a copy of the PMF of the interarrival
      /// time. 
      ///
      /// @return a copy of the pmf related to the interarrival time.
      PrositAux::pmf *getInterarrivalTime() const {
        return Z.get();
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

      /// @brief Check the periodicity of the task.
      ///
      /// This function returns a boolean with the periodicity of the task.
      ///
      /// @return the periodicity of the task.
      bool isPeriodic() const {
        return is_periodic;
      }

      /// @brief Obtain the task period.
      ///
      /// This function returns the period of the task. If the task is not
      /// periodic, an exception is raised.
      ///
      /// @return the task period.
      uint64_t getPeriod() const {
        return period;
      }

      /// @brief Obtain the deadline of the task.
      ///
      /// This function returns the deadline of the task.
      ///
      /// @return the deadline of the task.
      uint64_t getDeadline() const {
        return deadline;
      }

      /// @brief Set the deadline step.
      ///
      /// This function allows the user to decide the size of the deadline 
      /// step to be used in defining the probability map. All the deadlines 
      /// in the map are intended as multiples of this basic unit. The 
      /// operation generates an exception if the probabilistic deadlines 
      /// map is not empty.
      /// Additional constraints can be required by specific classes of 
      /// tasks, hence the need to make this function virtual.
      ///
      /// @param ds is the required deadline step.
      virtual void setDeadlineStep(DeadlineUnit ds);

      /// @brief Obtain the deadline step.
      ///
      /// This function allows the user to know the size of the deadline 
      /// step to be used in defining the probability map. All the deadlines 
      /// set in the map are intended as multiples of this basic unit.
      ///
      /// @return the current deadline step.
      DeadlineUnit getDeadlineStep() const {
        return deadline_step;
      }

      /// @brief Define a deadline in the map.
      ///
      /// This function inserts a deadline in the map of probabilistic 
      /// deadlines. An exception is thrown for duplicated entries and 
      /// when the specified deadline is not a multiple of the deadline
      /// step.
      ///
      /// @param the deadline to be inserted.
      void insertDeadline(DeadlineUnit deadline);

      /// @brief Obtain a pointer to the probabilistic deadlines map.
      ///
      /// This function returns a pointer to the probabilistic deadlines map.
      /// It is supposed to be used by the solver, hence the user is 
      /// discouraged to use this method directly.
      ///
      /// @return a pointer to the probabilistic deadlines map.
      DeadlineProbabilityMap *getProbabilisticDeadlines() {
        return &probabilistic_deadlines;
      }

      /// @brief Clear the deadline probability map.
      ///
      /// This function removes all previously set deadlines.
      void resetProbabilisticDeadlines() {
        probabilistic_deadlines.clear();
      }
      
      /// @brief Clear the deadline probability map.
      ///
      /// This function removes all previously set deadlines.
      void cleanProbabilisticDeadlines() {

        DeadlineProbabilityMapIter it = probabilistic_deadlines.begin();

        /* Iterate over the map */
        while (it != probabilistic_deadlines.end()) {

          /* Set the probability to zero */
          it->second = 0.0;

          /* Increment the ierator to point to next entry */
          it++;

        }

      }

/*      bool my_compare(const GenericTaskDescriptor &a, const GenericTaskDescriptor &b) {
        return a.deadline_step < b.deadline_step;
      }
*/

      virtual void setMinQoSTarget(double) { }
      virtual void setMaxQoSTarget(double) { }

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
          int index) = 0;

  };

}

#endif
