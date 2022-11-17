/*!
 * @file    help_messages.hpp
 * 
 * @brief   This file defines an implementation for a collection of 
 *          help messages presented in the code.
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
#ifndef HELP_MESSAGES_HPP
#define HELP_MESSAGES_HPP

#include "exceptions.hpp"

namespace PrositAux {

  /// @brief Prints the help for the CLI probability solver.
  ///
  /// This function is used to print the help function of the CLI solver.
  void help_cli();

  /// @brief Prints the help for the XML probability solver.
  ///
  /// This function is used to print the help function of the XML solver.
  void help_xml();

  /// @brief Prints the help for the autocorrelation program.
  ///
  /// This function is used to print the help function of the autocorrelation
  /// program.
  void help_autocorr();

  /// @brief Prints the help for the HMM Learning problem.
  ///
  /// This function is used to print the help function of the HMM Learner.
  void help_learner();

  /// @brief Prints the help for the HMM Decoding problem.
  ///
  /// This function is used to print the help function of the HMM Decoder.
  void help_decoder();

  /// @brief Prints the help for the CBS simulator.
  ///
  /// This function is used to print the help function of the CBS simulator.
  void help_simul();

}

#endif
