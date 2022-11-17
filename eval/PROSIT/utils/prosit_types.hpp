/*!
 * @file    prosit_types.hpp
 * 
 * @brief   Basic data types used in the PROSIT tool.
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
#ifndef PROSIT_TYPES_HPP
#define PROSIT_TYPES_HPP

#include <algorithm>     ///< std::min_element, std::min, std::max
#include <cassert>       ///< std::assert
#include <memory>        ///< std::unique_ptr
#include <vector>        ///< std::vector
#include <random>        ///< std::random_device
#include <cmath>         ///< std::abs
#include <list>          ///< std::list
#include <map>           ///< std::map

#include <unistd.h>

namespace PrositCore {

  /* Forward declaration */
  class FixedPriorityTaskDescriptor;

  /// @brief Base deadline.
  ///
  /// All deadlines are defined as integer multiples of a basic quantity.
  typedef uint64_t DeadlineUnit;

  /// @brief Deadline probability map.
  ///
  /// Data structure that associates a deadline with its probability.
  typedef std::map<DeadlineUnit, double> DeadlineProbabilityMap;

  /// @brief Iterator
  ///
  /// Iterator for the DeadlineProbabilityMap.
  typedef std::map<DeadlineUnit, double>::iterator DeadlineProbabilityMapIter;

  /// @brief Fixed priority tasks set.
  ///
  /// Data structure defining a set of fixed priority tasks.
  typedef std::vector<FixedPriorityTaskDescriptor*> FixedPriorityTaskSet;

  /// @brief Fixed priority tasks map.
  ///
  /// Data structure that associates a fixed priority with a set of tasks.
  typedef std::map<uint64_t, FixedPriorityTaskSet> FixedPriorityTaskMap;

  /// @brief Iterator
  ///
  /// Iterator for the FixedPriorityTasksMap.
  typedef std::map<uint64_t, FixedPriorityTaskSet>::iterator FixedPriorityTaskMapIter;

  /// @brief Structure for the results.
  struct results_t {

    int deadline;     ///< Relative deadline to be respected.
    int num_hits;     ///< Number of jobs that respect the deadline.
    int num_jobs;     ///< Total number of jobs.

    results_t() : deadline(0), num_hits(0), num_jobs(0) { }

  };

  /// @brief Structure for the response times.
  struct times_t {

    int64_t execution;     ///< Execution time of the job.
    uint64_t start;        ///< Starting time of the job.
    uint64_t end;          ///< Finishing time of the job.
    uint64_t response;     ///< Response time of the job.

    times_t() : execution(0), start(0), end(0), response(0) { }

  };

}

#endif
