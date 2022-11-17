/*!
 * @file    distribution.hpp
 * 
 * @brief   This class defines the header for probability distribution classes.
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
#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include <iomanip>         ///< std::setprecision
#include <fstream>         ///< std::ifstream
#include <Eigen/Dense>     ///< Eigen::

#include "exceptions.hpp"

namespace PrositAux {

  class distribution {

    protected:

      uint64_t min_val;           ///< The minimum index of the distribution.
      uint64_t max_val;           ///< The maximum index of the distribution.
      //uint64_t granularity;       ///< The granularity to represent the distribution.
      uint64_t size;              ///< The size of the distribution.
      Eigen::VectorXd elements;   ///< The values of the distribution.

    public:

      double epsilon;      ///< The epsilon error for the distribution.

      /// @brief Constructor.
      ///
      /// This is the constructor of the distribution.
      ///
      /// @param cmin is the miximum index of the distribution.
      /// @param cmax is the maximum index of the distribution.
      /// @param delta is the granularity of the distribution.
      /// @param eps is the epsilon error of the distribution.
      distribution(uint64_t cmin, uint64_t cmax/*, uint64_t delta*/, double eps) :
          min_val(cmin),
          max_val(cmax),
          //granularity(delta),
          size(cmax - cmin + 1),
          elements(cmax - cmin + 1),
          epsilon(eps) { }

      /// @brief Constructor.
      ///
      /// This is the copy constructor of the distribution.
      ///
      /// @param d is an already existing distribution.
      distribution(const distribution &d) :
          min_val(d.min_val),
          max_val(d.max_val),
          //granularity(d.granularity),
          size(d.size),
          elements(d.elements),
          epsilon(d.epsilon) { }

      /// @brief Destructor.
      ///
      /// This is the destructor of the distribution.
      virtual ~distribution() { }

      /// @brief Overloads the equality operator.
      ///
      /// This method overloads the equality operator for a distribution object.
      ///
      /// @param d is an already existing distribution.
      const distribution &operator=(const distribution &d) {

        min_val = d.min_val;
        max_val = d.max_val;
        //granularity = d.granularity;
        size = d.size;
        elements = d.elements;
        epsilon = d.epsilon;

        return *this;

      }

      /// @brief Return the maximum index of the distribution.
      ///
      /// This method returns the maximum element of the distribution.
      ///
      /// @return the maximum element of the distribution.
      uint64_t getMax() const {
        return max_val;
      }

      /// @brief Set the maximum index of the distribution.
      ///
      /// This method set the maximum element of the distribution.
      void setMax(uint64_t new_max) {
        max_val = new_max;
      }

      /// @brief Return the minimum index of the distribution.
      ///
      /// This method returns the minimum element of the distribution.
      ///
      /// @return the minimum element of the distribution.
      uint64_t getMin() const {
        return min_val;
      }

      /// @brief Set the minimum index of the distribution.
      ///
      /// This method sets the minimum element of the distribution.
      void setMin(uint64_t new_min) {
        min_val = new_min;
      }

      /// @brief Return the granularity of the distribution.
      ///
      /// This method returns the granularity of the distribution.
      ///
      /// @return the granularity of the distribution.
      /*uint64_t getGranularity() const {
        return granularity;
      }*/

      /// @brief Set the granularity of the distribution.
      ///
      /// This method sets the granularity of the distribution.
      /*void setGranularity(uint64_t new_granularity) {
        granularity = new_granularity;
      }*/

      /// @brief Return the size of the distribution.
      ///
      /// This method returns the number of element of the distribution.
      ///
      /// @return the size of the distribution.
      uint64_t getSize() const {
        return size;
      }

      /// @brief Set the size of the distribution.
      ///
      /// This method sets the number of element of the distribution.
      void setSize(uint64_t new_size) {
        size = new_size;
      }

      /// @brief Return the elements of the distribution.
      ///
      /// This method returns the elements of the distribution.
      ///
      /// @return the elements of the distribution.
      Eigen::VectorXd getElements() {
        return elements;
      }

      /// @brief Set the elements of the distribution.
      ///
      /// This method sets the elements of the distribution.
      void setElements(Eigen::VectorXd new_elements) {
        elements = new_elements;
      }

      /// @brief Load the file with the PMF of the distribution.
      ///
      /// This function is used to load the PMF of the distribution into
      /// the vector.
      ///
      /// @param filename is the file containing the distribution.
      ///
      /// @return zero if the load process was successful.
      virtual int load(const std::string &filename) = 0;

      /// @brief Print out the distribution.
      ///
      /// This function is used to print out the distribution.
      virtual void print() const;

      /// @brief Save the distribution.
      ///
      /// This function is used to save the distribution into a file.
      ///
      /// @param filename is the file containing the distribution.
      virtual void dump(const std::string &filename);

      /// @brief Set the value for an index in the distribution.
      ///
      /// This function is used to associate the value of the distribution
      /// with its respective index.
      ///
      /// @param idx is the index of the distribution.
      /// @param val is the value of the distribution.
      ///
      /// @return zero if the set process was successful.
      virtual int set(uint64_t idx, double val) = 0;

      /// @brief Get the value of the distribution.
      ///
      /// This function is used to obtain the value of the distribution
      /// associated to a given index.
      ///
      /// @param idx is the index of the distribution.
      ///
      /// @return the value of the distribution associated to the given index.
      virtual double get(uint64_t idx) const = 0;

  };

}

#endif
