/*!
 * @file    auxiliary_functions.hpp
 * 
 * @brief   This class defines the header for a collection of miscellaneous
 *          auxiliary mathematical functions.
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
#ifndef AUXILIARY_FUNCTIONS_HPP
#define AUXILIARY_FUNCTIONS_HPP

#include "prosit_types.hpp"
#include "exceptions.hpp"
#include "pmf.hpp"
#include "cdf.hpp"

namespace PrositAux {

  /// @brief Return the absolute time in microseconds.
  ///
  /// This function is used to read the system clock in order to take 
  /// perfromance metrics.
  long long my_get_time();

  /// @brief Check the size equality of two matrices.
  ///
  /// This function checks that two matrices are squared and that have equal 
  /// size. It is used to enforce that certain operations are legitimate.
  ///
  /// @param A0 is the first matrix to compare.
  /// @param A1 is the second matrix to compare.
  bool check_sizes(const Eigen::MatrixXd &A0, const Eigen::MatrixXd &A1);

  /// @brief Obtain the eigenvector related to the eigenvalue 1.
  ///
  /// This function returns the eigenvector related to eigenvalue 1.
  ///
  /// @param A is the matrix to analyze.
  Eigen::RowVectorXd stat(const Eigen::MatrixXd &A);

  /// @brief Load data into a matrix.
  ///
  /// This function loads the data from a file into a specified matrix.
  ///
  /// @param filename is the file containing the matrix.
  /// @param trans_mat is the transition matrix.
  void loadMatrix(const std::string &filename, Eigen::MatrixXd &trans_mat);

  /// @brief Return the maximum and minimum indices of a distribution.
  ///
  /// This function is used to determine the minimum and maximum indices
  /// of a distribution in order to create the proper vector.
  ///
  /// @param filename is the file containing the distribution.
  /// @param cmin is the minimum index of the distribution.
  /// @param cmax is the maximum index of the distribution.
  void obtainRange(const std::string &filename, uint64_t &cmin, uint64_t &cmax);

  /*! @brief Convert a PMF into a CDF */
  ///
  /// This function obtains the cumulative distribution function from a given
  /// probability mass function.
  ///
  /// @param prob is the PMF.
  /// @param cum is the CDF.
  void pmf2cdf(const pmf &prob, cdf &cum);

  /*! @brief Convert a CDF into a PMF */
  ///
  /// This function obtains the probability mass function from a given
  /// cumulative distribution function.
  ///
  /// @param cum is the CDF.
  /// @param prob is the PMF.
  void cdf2pmf(const cdf &cum, pmf &prob);

  /// @brief Obtain the infinity norm of a matrix.
  ///
  /// This function is used to compute the infinity norm of a matrix.
  ///
  /// @param A is the matrix to analyze.
  template <typename _Matrix_Type_>
  double InfinityNorm(const _Matrix_Type_ &A) {

    return A.cwiseAbs().rowwise().sum().maxCoeff();

  }

  /// @brief Obtain the Moore-Penrose pseudo-inverse of a matrix.
  ///
  /// This function is used to compute the Moore-Penrose pseudo-inverse of a 
  /// matrix.
  template <typename _Matrix_Type_>
  void pseudoInverse(const _Matrix_Type_ &a, _Matrix_Type_ &result,
      double epsilon = std::numeric_limits<double>::epsilon()) {

    Eigen::JacobiSVD<_Matrix_Type_> svd = 
        a.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);

    _Matrix_Type_ H(a.cols(), a.rows());

    H.setZero();
    double tolerance =
        double(epsilon * svd.singularValues().array().abs().maxCoeff());

    int len = svd.singularValues().array().rows();

    for (long i = 0; i < len; ++i) {

      double sigma = svd.singularValues()[i];

      if (sigma > tolerance)
        H(i, i) = typename _Matrix_Type_::Scalar(1.0) / 
            typename _Matrix_Type_::Scalar(sigma);

    }

    result = svd.matrixV() * H * (svd.matrixU().adjoint());

  }

  /// @brief Generate a random stochastic matrix.
  ///
  /// This function generates a random stochastic matrix.
  ///
  /// @param num_rows is the number of rows.
  /// @param num_cols is the number of columns.
  ///
  /// @return probability_matrix is the stochastic matrix.
  Eigen::MatrixXd generateProbabilityMatrix(int64_t num_rows, int64_t num_cols);

  /// @brief Load the observations samples.
  ///
  /// This function loads the file with the samples into a vector.
  ///
  /// @param filename is the name of the file to open.
  /// @param times is the vector with the samples.
  /// @param cmin is the minimum sample.
  /// @param cmax is the maximum sample.
  /// @param granularity is the granularity to resample the time.
  void loadObservations(const std::string filename, std::vector<uint64_t> &times, uint64_t &cmin, uint64_t &cmax, uint64_t granularity);

  /// @brief Load the emissions matrix.
  ///
  /// This function loads the file with the emissions matrix into a matrix.
  ///
  /// @param filename is the name of the file to open.
  /// @param emis_mat is the emissions matrix.
  /// @param cmin is the minimum observation sample.
  void loadEmissionMatrix(const std::string &filename, Eigen::MatrixXd &emis_mat, uint64_t cmin);

  /// @brief Save the probability matrix.
  ///
  /// This function saves the probability transition or emissions matrix into
  /// a file.
  ///
  /// @param filename is the name of the file to save.
  /// @param matrix is the matrix to be saved.
  /// @param cmin is the minimum observation sample.
  /// @param cmin is the maximum observation sample.
  void saveMatrix(const std::string &filename, const Eigen::MatrixXd &matrix, uint64_t cmin, uint64_t cmax);

  /// @brief Save the observations samples.
  ///
  /// This function saves the samples into a file.
  ///
  /// @param filename is the name of the file to save.
  /// @param times is the vector with the samples.
  void saveObservations(const std::string &filename, const std::vector<uint64_t> &times);

}

#endif
