/*!
 * @file    auxiliary_functions.cpp
 * 
 * @brief   This class defines an implementation for a collection of 
 *          miscellaneous auxiliary mathematical functions.
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
#include "auxiliary_functions.hpp"

using namespace std;
using namespace Eigen;

namespace PrositAux {

  /// @brief Return the absolute time in microseconds.
  long long my_get_time() {

    struct timespec time;

    /* Obtain the time */
    clock_gettime(CLOCK_MONOTONIC, &time);

    /* Convert it in microseconds */
    return (long long)((time.tv_sec * 1000000LL) + (time.tv_nsec / 1000LL));

  }

  /// @brief Check the size equality of two matrices.
  bool check_sizes(const MatrixXd &A0, const MatrixXd &A1) {

    return (A0.rows() == A0.cols()) && (A0.rows() == A1.rows()) &&
        (A0.cols() == A1.cols());

  }

  /// @brief Obtain the eigenvector related to the eigenvalue 1.
  RowVectorXd stat(const MatrixXd &A) {

    if (A.rows() != A.cols())
      throw Exc("Function Stat: The matrix A has to be square");

    VectorXd e(A.rows());
    MatrixXd B(A.rows(), A.cols() + 1);
    MatrixXd P;
    e.setOnes();

    B.block(0, 0, A.rows(), A.cols()) = A - 
        MatrixXd::Identity(A.rows(), A.cols());

    B.block(0, A.cols(), e.rows(), e.cols()) = e;

    RowVectorXd y(A.rows() + 1);
    y.setZero();
    y(A.rows()) = 1;
    pseudoInverse<MatrixXd>(B, P);

    return y * P;

  }

  /// @brief Load the transition matrix.
  void loadMatrix(const string &filename, MatrixXd &trans_mat) {

    /* Variable declaration */
    ifstream data;
    int64_t num_rows = trans_mat.rows();
    int64_t num_cols = trans_mat.cols();
  
    /* Open the data file */
    data.open(filename.c_str());
    if (data.fail()) {
      EXC_PRINT_2("ERROR: Unable to find the file: ", filename);
    }

    /* Iterate over the rows of the file */
    for (int64_t i = 0; i < num_rows; i++) {

      /* Iterate over the columns of the file */
      for (int64_t j = 0; j < num_cols; j++) {

        /* Fill the matrix */
        data >> trans_mat(i, j);

      }

    }

    /* Close the file */
    data.close();

  }

  /// @brief Load the emission matrix.
  void loadEmissionMatrix(const string &filename, MatrixXd &emis_mat, uint64_t cmin) {

    int64_t num_states = emis_mat.rows();

    /* Iterate over the number of modes */
    for (int64_t i = 0; i < num_states; i++) {

      /* Variables for the files */
      ifstream myfile;

      /* Variables for the data */
      int index;
      long double value;

      /* Variable for the name of the computation times */
      char comp_time_file[250];

      /* Build the string with the name */
      sprintf(comp_time_file, "%s%lu.txt", filename.c_str(), i + 1);

      /* Open the data file */
      myfile.open(comp_time_file);
      if (myfile.fail()) {
        EXC_PRINT_2("ERROR: Unable to find the file: ", comp_time_file);
      }

      /* Check if the file was properly open */
      if (!myfile.is_open()) {
        EXC_PRINT_2("ERROR: Unable to open the file: ", comp_time_file);
      }

      /* Iterate over the file */
      while (myfile >> index >> value) {

        /* Set the element in the matrix */
        emis_mat(i, index - cmin) = value;
  
      }

      /* Close the file */
      myfile.close();

    }

  }

  /// @brief Determine the range of the given distribution.
  void obtainRange(const string &filename, uint64_t &cmin, uint64_t &cmax) {

    /* Variables for the files */
    ifstream myfile;
    string line;

    /* Counter of elements within each line */
    int size;

    /* Variables for the data */
    double value;
    uint64_t index;
  
    /* Open the data file */
    myfile.open(filename.c_str());
    if (myfile.fail()) {
      EXC_PRINT_2("ERROR: Unable to find the file: ", filename);
    }

    /* Check if the file was properly open */
    if (!myfile.is_open()) {
      EXC_PRINT_2("ERROR: Unable to open the file: ", filename);
    }

    /* Iterate over the file */
    while (myfile.good()) {
 
      /* Obtain a line of the file */
      getline(myfile, line);
 
      /* Initialize the elements counter */
      size = 0;

      /* Check whether the lines are over */
      if (line.find_first_not_of(' ') == std::string::npos) {
        break;
      }

      /* Obtain the line */
      istringstream sstr1(line);

      /* Iterate over the line */
      while (!sstr1.eof()) {

        /* Check whether the element is numeric */
        if (!(sstr1 >> value)) {
          EXC_PRINT_2("ERROR: Incorrect format for file: ", filename);
        }

        /* Update the number of elements within the line */
        size++;

      }

      /* Check the correct number of elements within each line */
      if (size > 2) {
        EXC_PRINT_2("ERROR: Unknown format for file: ", filename);
      }

      /* Re-obtain the line */
      istringstream sstr2(line);

      /* Parse each element of the line */
      sstr2 >> index >> value;

      /* Check if the PMF has negative values */
      if (value < 0) {
        EXC_PRINT_2("ERROR: There are negative values in the file: ", filename);
	  }

      /* Determine the minimum index */
      if (index < cmin) {
        cmin = index;
      }

      /* Determine the maximum index */
      if (index > cmax) {
        cmax = index;
      }

    }

    /* Close the file */
    myfile.close();

  }

  /// @brief Convert a PMF into a CDF.
  void pmf2cdf(const pmf &prob, cdf &cum) {

    /* Initialize the variables */
    double sum = 0.0;

    /* Iterate over the PMF */
    for (uint64_t i = 0; i < prob.getSize(); i++) {

      /* Accumulate the probability */
      sum += prob.get(i);

      /* Set the value */
      cum.set(i, sum);

    }

    /* Complete the CDF to 1 */
    if (cum.get(cum.getMax() - cum.getMin()) < (1.0 - cum.epsilon)) {
      cerr << "WARNING: The resulting CDF did not sum to 1." << endl;
      cum.set(cum.getMax() - cum.getMin(), 1.0);
    }

    /* Check possible errors with the CDF */
    if (cum.check() != cdf::ERR_CODES::CDF_OK)
      EXC_PRINT("ERROR: The resulting CDF was ill formed.");

    return;

  }

  /// @brief Convert a CDF into a PMF.
  void cdf2pmf(const cdf &cum, pmf &prob) {

    /* Initialize the variables */
    double prev = 0.0;

    /* Iterate over the CDF */
    for (uint64_t i = 0; i < cum.getSize(); i++) {

      /* Set the value */
      prob.set(i, cum.get(i) - prev);

      /* Update the previous value */
      prev = cum.get(i);

    }

    return;

  }

  /// @brief Generate a random stochastic matrix.
  MatrixXd generateProbabilityMatrix(int64_t num_rows, int64_t num_cols) {

    /* Generate a positive valued random matrix */
    MatrixXd probability_matrix = (MatrixXd::Random(num_rows, num_cols)).cwiseAbs();

    /* Obtain the vector with the row sum */
    VectorXd sum_vector = (probability_matrix.rowwise().sum()).cwiseInverse();

    /* Normalize the probability matrix */
    probability_matrix = probability_matrix.array().colwise() * sum_vector.array();

    return probability_matrix;

  }

  /// @brief Load the observations samples.
  void loadObservations(const string filename, vector<uint64_t> &times, uint64_t &cmin, uint64_t &cmax, uint64_t granularity) {

    /* Variable declaration */
    ifstream data;
    uint64_t comp_time;
  
    /* Open the data file */
    data.open(filename.c_str());
    if (data.fail()) {
      EXC_PRINT_2("ERROR: Unable to find the file: ", filename);
    }

    /* Iterate over the file */
    while (data >> comp_time) {

      /* Resample the time to the desired granularity */
      comp_time = (comp_time + granularity - 1) / granularity;

      /* Fill the vector */
      times.push_back(comp_time);

      /* Determine the minimum value */
      if (comp_time < cmin) {
        cmin = comp_time;
      }

      /* Determine the maximum value */
      if (comp_time > cmax) {
        cmax = comp_time;
      }

    }
    //cout << cmin << " " << cmax << endl;

    /* Close the file */
    data.close();

  }

  /// @brief Save the matrix.
  void saveMatrix(const string &filename, const MatrixXd &matrix, uint64_t cmin, uint64_t cmax) {

    /* Variables for the files */
    ofstream dumpFile;

    /* Open the data file */
    dumpFile.open(filename.c_str());

    /* Check if the file was properly open */
    if(dumpFile.is_open()) {

      if ((cmin == 0) && (cmax == 0)) {

        /* Iterate over the rows */
        for (int64_t i = 0; i < matrix.rows(); i++) {

          /* Iterate over the columns */
          for (int64_t j = 0; j < matrix.cols(); j++) {

            /* Save the matrix */
            dumpFile << fixed << setprecision(14) << matrix(i, j) << " ";

          }

          /* Jump to the next row */
          dumpFile << endl;

        }

      }

      else {

        /* Iterate over the rows */
        for (uint64_t i = cmin; i <= cmax; i++) {

          /* Save the vector */
          dumpFile << fixed << i << " " << setprecision(14) << matrix(i - cmin) << endl;

        }

      }

      /* Close the data file */
      dumpFile.close();

    }
    else {

      cerr << "WARNING: Unable to open the file: " << filename << endl;

    }

  }

  /// @brief Save the observations samples.
  void saveObservations(const string &filename, const vector<uint64_t> &times) {

    /* Variables for the files */
    ofstream dumpFile;

    /* Open the data file */
    dumpFile.open(filename.c_str());

    /* Check if the file was properly open */
    if(dumpFile.is_open()) {

      /* Iterate over the vector of times */
      for (size_t i = 0; i < times.size(); i++) {

        /* Prints out the current element */
        dumpFile << times[i] << endl;

      }

      /* Close the data file */
      dumpFile.close();

    }
    else {

      cerr << "WARNING: Unable to open the file: " << filename << endl;

    }

  }

}
