/*!
 * @file    cdf.hpp
 * 
 * @brief   This class defines the header for cumulative distribution function
 *          classes.
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
#ifndef CDF_HPP
#define CDF_HPP

#include "distribution.hpp"

namespace PrositAux {

  /* Forward declaration */
  class pmf;

  class cdf : public distribution {

    private:

      bool just_created;     ///< The .

    public:

      /// @brief Structure for the distribution error codes.
      ///
      /// This structure contains the different error codes that can occur
      /// when loading the CDF.
      typedef enum ERR_CODES {
        CDF_OK,
        CDF_BAD_MAX,
        CDF_NON_INCREASING
      } ERR_CODES;

      /// @brief Constructor.
      ///
      /// This is the constructor of the distribution.
      ///
      /// @param cmin is the miximum index of the distribution.
      /// @param cmax is the maximum index of the distribution.
      /// @param delta is the granularity of the distribution.
      /// @param eps is the epsilon error of the distribution.
      cdf(int cmin, int cmax/*, int delta*/, double eps);

      /// @brief Constructor.
      ///
      /// This is the copy constructor of the distribution.
      ///
      /// @param d is an already existing distribution.
      cdf(const cdf &d) : distribution::distribution(d) {
        just_created = d.just_created;
      }

      /// @brief Destructor.
      ///
      /// This is the destructor of the distribution.
      virtual ~cdf() { }

      /// @brief Check the CDF of the distribution.
      ///
      /// This function is used to check the CDF of the distribution to 
      /// determine whether is well-formed.
      ///
      /// @return the error code obtained from the check process.
      ERR_CODES check() const;

      /// @brief Load the file with the CDF of the distribution.
      ///
      /// This function is used to load the CDF of the distribution into
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
