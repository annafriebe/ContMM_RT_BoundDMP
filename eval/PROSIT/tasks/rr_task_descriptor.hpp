/*!
 * @file    rr_task_descriptor.hpp
 * 
 * @brief   This class defines the header for the resource reservation task 
 *          descriptor, where the tasks are scheduled using resource 
 *          reservations.
 *          An obvious choice for resource reservations is to have deadlines 
 *          that are integer multiples of the server period.
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
#ifndef RR_TASK_DESCRIPTOR_HPP
#define RR_TASK_DESCRIPTOR_HPP

#include "../solver/rr_probability_solver.hpp"
#include "generic_task_descriptor.hpp"

namespace PrositCore {

  class ResourceReservationTaskDescriptor : public GenericTaskDescriptor {

    protected:

      bool solved;        ///< The solver has been called.
      int num_modes;      ///< Number of modes.

      ///< Distribution of the computation time.
      std::vector<std::unique_ptr<PrositAux::pmf>> C;

      ///< Probability transition matrix for the mode changes.
      Eigen::MatrixXd transition_matrix;

      uint64_t Qs;             ///< Reservation budget.
      uint64_t Ts;             ///< Reservation period.
      uint64_t granularity;    ///< Granularity for the resampling.

      ///< Pointer to the object containing the solution algorithm for 
      ///< probabilities.
      std::unique_ptr<ResourceReservationProbabilitySolver> probability_solver;

      uint64_t Qmin;
      uint64_t Qmax;
      double QoSmin;            ///< Reservation budget.
      double QoSmax;            ///< Reservation budget.
      double qos_min_target;    ///< Reservation budget.
      double qos_max_target;
      bool bounds_inited;

    public:

      double bandwidth = 0.0;    ///< Total bandwidth used.

      /// @brief Constructor.
      ///
      /// This is the constructor for resource reservation tasks.
      ///
      /// @param nm is the unique identifier for the task.
      /// @param num_modesd is the number of modes for the task.
      /// @param Cd is the vector of pointers to the distribution of the computation time.
      /// @param Zd is the pointer to the distribution of the interarrival time.
      /// @param transition_matrixd is the transition probability matrix for the mode change.
      /// @param Qd is the reservation budget.
      /// @param Tsd is the reservation period.
      /// @param algorithm is the solver to be set.
      ResourceReservationTaskDescriptor(std::string nm, 
          int num_modesd,
          std::vector<std::unique_ptr<PrositAux::pmf>> Cd, 
          std::unique_ptr<PrositAux::pmf> Zd, 
          Eigen::MatrixXd transition_matrixd,
          const uint64_t Qsd, 
          const uint64_t Tsd,
          uint64_t granularityd, 
          uint64_t deadlined,
          std::string algorithm);

      /// @brief Destructor.
      ///
      /// This is the destructor of the resource reservation task descriptor.
      virtual ~ResourceReservationTaskDescriptor();

      /// @brief Obtain the number of modes of the task.
      ///
      /// This function returns the number of modes of the task.
      ///
      /// @return the number of modes of the task.
      int getNumModes() const {
        return num_modes;
      }

      /// @brief Obtain a copy of the computation time distribution.
      ///
      /// This function is used to return a copy of the PMF of the computation 
      /// time.
      ///
      /// @return a copy of the pmf related to the computation time.
      std::vector<std::unique_ptr<PrositAux::pmf>> const& getComputationTime() const {
        return C;
      }

      /// @brief Obtain a copy of the transition matrix.
      ///
      /// This function is used to return a copy of the transition matrix 
      /// of the mode change.
      ///
      /// @return a copy of the transition matrix associated to the task.
      Eigen::MatrixXd getTransitionMatrix() {
        return transition_matrix;
      }

      /// @brief Sets the reservation budget.
      ///
      /// This function sets the budget of the task. An exception is thrown if 
      /// the bandwidth exceeds 1.0
      ///
      /// @param Qd is the desired budget. 
      void setBudget(uint64_t Qd) {

        /* Check if the bandwidth is smaller than 1 */
        if ((double(Qd) / double(Ts)) > 1.0) {
          EXC_PRINT_2("ERROR: The budget is too large for task ", name);
        }

        /* Set the budget */
        Qs = Qd;

        /* Update the bandwidth */
        bandwidth = double(Qs) / double(Ts);

      }

      /// @brief Obtain the reservation budget.
      ///
      /// This function returns the budget of the task.
      ///
      /// @return the current budget.
      uint64_t getBudget() const {
        return Qs;
      }

      /// @brief Sets the server period
      ///
      /// @param Tsd desired period. AN exception is thrown if the
      /// bandwidth exceeds 1.0 or if the server period is not a sub multiple
      /// of current deadline step.
      void setServerPeriod(uint64_t Tsd) {

        
        /* Check if the bandwidth is smaller than 1 */
        if ((double(Qs) / double(Tsd)) > 1.0) {
          EXC_PRINT_2("ERROR: The server period is too small for task ", name);
        }

        /* Check if the deadline step is a multiple of the server period */
        if ((deadline_step != 0) && (!(deadline_step % Tsd))) {
          EXC_PRINT_2("ERROR: The deadline step has to be a multiple of the "
              "server period for ", name);
        }

        /* Set the reservation period */
        Ts = Tsd;

        /* Update the bandwidth */
        bandwidth = double(Qs) / double(Ts);

      }

      /// @brief Obtain the reservation period.
      ///
      /// This function returns the reservation period of the task.
      ///
      /// @return the current reservation period.
      uint64_t getServerPeriod() const {
        return Ts;
      }

      void setGranularity(uint64_t granularityd) {
        granularity = granularityd;
      }

      uint64_t getGranularity() const {
        return granularity;
      }

      /// @brief Set the deadline step.
      ///
      /// This function overrides the standard definition of deadline step,
      /// requiring the additional constraint that it has to be a multiple 
      /// of the server period. It can be freely reset to zero at any time.
      ///
      /// @param ds is the required deadline step.
      void setDeadlineStep(DeadlineUnit ds);

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
      void setSolver(std::unique_ptr<ResourceReservationProbabilitySolver> psd);

      void resetSolver();

      /// @brief Compute the probability of respecting the deadlines.
      ///
      /// These probabilities are computed for a given configuration of the
      /// scheduling parameters. The deadlines have been previously registered
      /// by calling the function insertDeadline(). An exception is thrown if 
      /// the solver is not set.
      void computeProbability();

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
      double getProbability(DeadlineUnit deadline);

      /// @brief Evaluate the QoS function for a given budget.
      double QoSEvaluation(DeadlineUnit Q);

      bool inverseQoSEvaluation(double p, uint64_t &Q, bool ceil);

      bool identifyBounds(double q_min, double q_max);

      uint64_t getQMin() const {
        return Qmin;
      }

      uint64_t getQMax() const {
        return Qmax;
      }

      uint64_t setQMin(int Qmind) {

        uint64_t QMin_old = Qmin;
        Qmin = Qmind;
        bounds_inited = false;
        return QMin_old;

      }

      uint64_t setQMax(int Qmaxd) {

        uint64_t QMax_old = Qmaxd;
        Qmax = Qmaxd;
        bounds_inited = false;
        return QMax_old;

      }

      double getQoSMin() const {
        return QoSmin;
      }

      double getQoSMax() const {
        return QoSmax;
      }

      double getMinQoSTarget() const {
        return qos_min_target;
      }

      double getMaxQoSTarget() const {
        return qos_max_target;
      }

      void setMinQoSTarget(double q) {

        qos_min_target = q;
        bounds_inited = false;

      }

      void setMaxQoSTarget(double q) {

        qos_max_target = q;
        bounds_inited = false;

      }

      bool getBoundsInited() {
        return bounds_inited;
      }

      void initBounds() {

        QoSmin = QoSEvaluation(Qmin);
        QoSmax = QoSEvaluation(Qmax);
        bounds_inited = true;

      }

      /// @brief Displays results for a single task.
      ///
      /// @param td is the ...
      /// @param probability is the ...
      /// @param quality is the ...
      /// @param time is the ...
      /// @param show_time is the ...
      /// @param index is the ...
      virtual void display(std::string,
          GenericTaskDescriptor* td, 
          const std::vector<double> &probability, 
          const std::vector<double> &quality, 
          const std::vector<long long> &time, 
          int index);

  };

}

#endif
