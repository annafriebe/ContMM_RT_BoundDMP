/*!
 * @file    pmf.hpp
 * 
 * @brief   This class defines the header for probability mass function classes.
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
#ifndef PMF_HPP
#define PMF_HPP

#include <memory>     ///< std::unique_ptr

#include "distribution.hpp"

namespace PrositAux {

  /* Forward declaration */
  class cdf;

  class pmf : public distribution {

    private:

      double tail;     ///< The tail of the distribution.

    public:

      /// @brief Structure for the distribution error codes.
      ///
      /// This structure contains the different error codes that can occur
      /// when loading the PMF.
      typedef enum ERR_CODES {
        PMF_OK,
        PMF_SMALLER_ONE,
        PMF_GREATER_ONE
      } ERR_CODES;

      /// @brief Constructor.
      ///
      /// This is the constructor of the distribution.
      ///
      /// @param cmin is the miximum index of the distribution.
      /// @param cmax is the maximum index of the distribution.
      /// @param delta is the granularity of the distribution.
      /// @param eps is the epsilon error of the distribution.
      pmf(uint64_t cmin, uint64_t cmax/*, uint64_t delta*/, double eps) :
          distribution(cmin, cmax/*, delta*/, eps), tail(0.0) {

        /* Initialize the vector with the distribution */
        elements.setZero();

      }

      /// @brief Constructor.
      ///
      /// This is the copy constructor of the distribution.
      ///
      /// @param d is an already existing distribution.
      pmf(const pmf &d) : distribution::distribution(d) {
        tail = d.tail; 
      }

      /// @brief Destructor.
      ///
      /// This is the destructor of the distribution.
      virtual ~pmf() { }

      /// @brief Return the tail of the distribution.
      ///
      /// This method returns the tail of the distribution.
      ///
      /// @return the tail of the distribution.
      double getTail() {
        return tail;
      }

      /// @brief Obtain the sum of the PMF.
      ///
      /// This function is used to obtain the sum of the distribution.
      ///
      /// @return the sum of the distribution.
      double sum() const;

      /// @brief Obtain the mean of the PMF.
      ///
      /// This function is used to obtain the mean of the distribution.
      ///
      /// @return the mean of the distribution.
      double avg() const;

      /// @brief Obtain the standard deviation of the PMF.
      ///
      /// This function is used to obtain the standard deviation of the 
      /// distribution.
      ///
      /// @return the standard deviation of the distribution.
      double std() const;

      /// @brief Obtain the variance of the PMF.
      ///
      /// This function is used to obtain the variance of the distribution.
      ///
      /// @return the variance of the distribution.
      double var() const;

      /// @brief Check the PMF of the distribution.
      ///
      /// This function is used to check the PMF of the distribution to 
      /// determine whether is well-formed.
      ///
      /// @return the error code obtained from the check process.
      ERR_CODES check() const;

      /// @brief Load the file with the PMF of the distribution.
      ///
      /// This function is used to load the PMF of the distribution into
      /// the vector.
      ///
      /// @param filename is the file containing the distribution.
      ///
      /// @return zero if the load process was successful.
      virtual int load(const std::string &filename);

      /// @brief Set the value for an index in the distribution.
      ///
      /// This function is used to associate the value of the distribution
      /// with its respective index.
      ///
      /// @param idx is the index of the distribution.
      /// @param val is the value of the distribution.
      ///
      /// @return zero if the set process was successful.
      virtual int set(uint64_t idx, double val);

      /// @brief Get the value of the distribution.
      ///
      /// This function is used to obtain the value of the distribution
      /// associated to a given index.
      ///
      /// @param idx is the index of the distribution.
      ///
      /// @return the value of the distribution associated to the given index.
      virtual double get(uint64_t idx) const;

      /// @brief Resample the PMF to a given granularity
      ///
      /// Creates a clone of the PMF, which is resampled at the requested
      /// granularity, always being conservative. This means that resampling
      /// the computation time must aproximate to the greater value (ceil),
      /// while with the interarrival time the aproximation must be done
      /// towards a lower value (floor). This strategy ensures an always
      /// conservative approach.
      ///
      /// This is what the growing flag is used for:
      ///   True  -> Computation  time -> greater value.
      ///   False -> Interarrival time -> lower value.
      ///
      /// @param delta is the granularity for the resampling.
      /// @param grow is the direction of the aproximation (left/right).
      ///
      /// @return the resampled clone of the distribution.
      std::unique_ptr<pmf> resample(uint64_t delta, bool growing);

      /// @brief Convert a PMF into a CDF.
      ///
      /// This function is used to transform a given PMF into a CDF.
      ///
      /// @param prob is the PMF.
      /// @param cum is the CDF.
      friend void pmf2cdf(const pmf &prob, cdf &cum);

      /// @brief Convert a CDF into a PMF.
      ///
      /// This function is used to transform a given CDF into a PMF.
      ///
      /// @param cum is the CDF.
      /// @param prob is the PMF.
      friend void cdf2pmf(const cdf &cum, pmf &prob);

  };

}

#endif
