/*!
 * @file    qbd_rr_probability_solver.hpp
 * 
 * @brief   This class defines the header for a family of solvers for resource 
 *          reservation tasks based on algorithms for numeric solution of 
 *          Quasi-Birth-Death processes.
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
#ifndef QBD_RR_PROBABILITY_SOLVER_HPP
#define QBD_RR_PROBABILITY_SOLVER_HPP

#include "../utils/auxiliary_functions.hpp"
#include "rr_probability_solver.hpp"

namespace PrositCore {

  class QBDResourceReservationProbabilitySolver : public ResourceReservationProbabilitySolver {

    protected:

      Eigen::MatrixXd B0;      ///< Upper-left submatrix of the QBDP transition matrix.
      Eigen::MatrixXd A0;      ///< Lower-left submatrix of the QBDP transition matrix.
      Eigen::MatrixXd A1;      ///< Lower-right submatrix of the QBDP transition matrix.
      Eigen::MatrixXd A2;      ///< Upper-right submatrix of the QBDP transition matrix.
      Eigen::MatrixXd R;       ///< Minimal nonnegative solution to R = A2 + R A1 + R^2 A0.
      Eigen::RowVectorXd pi0;  ///< Stationary vector of the QBDP.

      /// @brief Prepare the computation.
      ///
      /// This function generates the matrices used by the QBD methods.
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

      /// @brief Computes the stationary vector pi0.
      ///
      /// This function computes the stationary vector pi0 of the QBDP from the
      /// obtained matrices R and A0.
      bool computePi0();

    public:

      /// @brief Constructor.
      ///
      /// This is the constructor for QBD resource reservation probability 
      /// solvers.
      QBDResourceReservationProbabilitySolver() : 
          ResourceReservationProbabilitySolver(),
          B0(),
          A0(),
          A1(),
          A2(),
          R(),
          pi0() { }

      /// @brief Destructor.
      ///
      /// This is the destructor for QBD resource reservation probability 
      /// solvers.
      virtual ~QBDResourceReservationProbabilitySolver() { }

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
