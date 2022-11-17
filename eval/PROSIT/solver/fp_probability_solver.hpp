/*!
 * @file    fp_probability_solver.hpp
 * 
 * @brief   This class contains the abstract definitions for generic 
 *          probability solvers for a set of fixed priority tasks. In 
 *          essence, a probability solver wraps the different algorithms 
 *          used to compute probabilistic deadlines. The solver for fixed 
 *          priority tasks considers sets of tasks.
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
#ifndef FP_PROBABILITY_SOLVER_HPP
#define FP_PROBABILITY_SOLVER_HPP

#include "../utils/auxiliary_functions.hpp"
#include "../tasks/fp_task_schedule.hpp"
#include "probability_solver.hpp"

namespace PrositCore {

  /* Forward declaration */
  class FixedPriorityTaskSchedule;

  class FixedPriorityProbabilitySolver : public ProbabilitySolver {

    protected:

      double epsilon;    ///< Threshold for the result of two iterations.
      int max_iter;      ///< Maximum number of iteration of the algorithm.

      Eigen::MatrixXd B0;      ///< Upper-left submatrix of the QBDP transition matrix.
      Eigen::MatrixXd A0;      ///< Lower-left submatrix of the QBDP transition matrix.
      Eigen::MatrixXd A1;      ///< Lower-right submatrix of the QBDP transition matrix.
      Eigen::MatrixXd A2;      ///< Upper-right submatrix of the QBDP transition matrix.
      Eigen::MatrixXd R;       ///< Minimal nonnegative solution to R = A2 + R A1 + R^2 A0.
      Eigen::RowVectorXd pi0;  ///< Stationary vector of the QBDP.

      ///< The task descriptors are ordered by their scheduling priority.
      FixedPriorityTaskSchedule *task_schedule;

      std::unique_ptr<PrositAux::pmf> addDistribution(std::unique_ptr<PrositAux::pmf> &first, const std::unique_ptr<PrositAux::pmf> &second);

      void normalizeDistribution(std::unique_ptr<PrositAux::pmf> &dist);

      void shiftDistribution(std::unique_ptr<PrositAux::pmf> &dist, uint64_t delta);

      std::unique_ptr<PrositAux::pmf> computeBacklogUntil(const std::unique_ptr<PrositAux::pmf> &starting_distribution, uint64_t id_task, uint64_t until_time);

      std::unique_ptr<PrositAux::pmf> splitConvolve(std::unique_ptr<PrositAux::pmf> &first, const std::unique_ptr<PrositAux::pmf> &second, uint64_t delta_second);

      std::unique_ptr<PrositAux::pmf> convolveShrink(std::unique_ptr<PrositAux::pmf> &first, const std::unique_ptr<PrositAux::pmf> &second, uint64_t delta_second);

      std::unique_ptr<PrositAux::pmf> convolutionsShrink(std::unique_ptr<PrositAux::pmf> &dist, uint64_t delta);

      /// @brief Computes the response time for a job.
      std::unique_ptr<PrositAux::pmf> getJobResponseTime(std::unique_ptr<PrositAux::pmf> &job_back, uint64_t id_task, uint64_t curr_activation);

      /// @brief Perform the convolution of two given distributions.
      ///
      /// This function performs the convolution of two given distributions.
      ///
      /// @param first is the input distribution.
      /// @param second is the impulse response distribution.
      ///
      /// @return clone is the convolved distribution.
      std::unique_ptr<PrositAux::pmf> convolutionsConvolve(const std::unique_ptr<PrositAux::pmf> &first, const std::unique_ptr<PrositAux::pmf> &second);

      /// @brief Check the conditions to start.
      ///
      /// This function determines whether everything is ok to start the 
      /// pre-process phase.
      bool checkList();

      /// @brief Prepare the computation.
      ///
      /// This function generates the matrices used by the QBDP methods.
      void preProcess();

      /// @brief Post processing after the solution.
      ///
      /// This function triggers the computation of the stationary vector of
      /// the QBDP.
      void postProcess();

      /// @brief Fills in the probability map.
      ///
      /// This function computes the probability map from the obtained 
      /// stationary vector pi0.
      void fillInProbabilityMap();

      /// @brief Executes the iterations of the Logarithmic Reduction method.
      ///
      /// This function implements the Logarithmic Reduction method for discrete
      /// time QBDP. The iteration produces a matrix R that has to be 
      /// post-processed.
      void applyAlgorithm();

      /// @brief Computes the stationary vector pi0.
      ///
      /// This function computes the stationary vector pi0 of the QBDP from the
      /// obtained matrices R and A0.
      bool computePi0();

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for fixed priority probability solvers.
      FixedPriorityProbabilitySolver(double epsilon_d, unsigned int max_iter_d) : 
          ProbabilitySolver(), 
          epsilon(epsilon_d), 
          max_iter(max_iter_d),
          task_schedule(0) {

        /* Check the correctness of the epsilon */
        if (epsilon < 0) {
          EXC_PRINT("ERROR: The epsilon parameter has to be non negative.");
        }

      }

      /// @brief Destructor.
      ///
      /// This is the destructor for fixed priority probability solvers.
      virtual ~FixedPriorityProbabilitySolver() { }

      /// @brief Associates the solver with a task schedule.
      ///
      /// A fixed priority solver applies to a fixed priority schedule.
      ///
      /// @param ts is the descriptor of the schedule.
      virtual void registerSchedule(FixedPriorityTaskSchedule *ts);

    private:

      /// @brief Extract the QBDP submatrices.
      ///
      /// This function extracts the submatrices from a QBDP having a recursive
      /// structure. The QBPD has the form:
      ///
      ///          B  A2  0  0  0 ...
      ///          A0 A1 A2  0  0 ...
      ///     P =  0  A0 A1 A2  0 ...
      ///          0  0  A0 A1 A2 ...
      ///          ...
      ///
      /// @param mat is the probability transition matrix of the QBDP.
      /// @param dim is the size of the submatrices.
      /// @param B  is the submatrix B.
      /// @param A0 is the submatrix A0.
      /// @param A1 is the submatrix A1.
      /// @param A2 is the submatrix A2.
      void extractSubmatrices(const Eigen::MatrixXd &P,
          uint64_t dim,
          Eigen::MatrixXd &B,
          Eigen::MatrixXd &A0,
          Eigen::MatrixXd &A1,
          Eigen::MatrixXd &A2);

  };

}

#endif
