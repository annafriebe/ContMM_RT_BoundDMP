/*!
 * @file    exceptions.hpp
 * 
 * @brief   This class defines the exception objects used throughout the library.
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
#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

#include <iostream>
#include <stdexcept>

namespace PrositAux {

  class Exc: public std::runtime_error {

    public:

      /// @brief Constructor.
      ///
      /// @param msg is the error message.
      Exc(const char *ms) : std::runtime_error(ms) { }

      Exc(const std::string &s) : std::runtime_error(s) { }

  };




  /// @brief Macro to print out an error message and throw an exception at the same 
  ///        time.
  ///
  /// If NDEBUG is defined, the macro simply throws the exception. Otherwise, it 
  /// prints out the msg, along with the file, line and function where the exception
  /// occurred.
  ///
  /// @param msg is the error message.
  #ifdef NDEBUG
    #define EXC_PRINT(msg) throw Exc(std::string(__PRETTY_FUNCTION__)+std::string(": ")+std::string(msg))
  #else
    #define EXC_PRINT(msg)                                                                                          \
    do {                                                                                                            \
      std::cerr<<"file:"<<__FILE__<<", line:"<<__LINE__<<", function: "<<__PRETTY_FUNCTION__<<". "<<std::endl<<msg<<std::endl; \
      throw PrositAux::Exc(std::string(__PRETTY_FUNCTION__)+std::string(": ")+std::string(msg)); \
    } while(0)
  #endif




  /// @brief Macro to print out two error messages and throw an exception (with the
  ///        message stored) at the same time.
  ///
  /// If NDEBUG is defined, the macro simply throws the exception. Otherwise, it
  /// prints out the two messages along with the file, line and function where
  /// the exception occurred.
  ///
  /// @param msg1 is the first error message.
  /// @param msg2 is the second error message.
  #ifdef NDEBUG 
    #define EXC_PRINT_2(msg1, msg2) throw Exc(std::string(__PRETTY_FUNCTION__)+std::string(": ")+std::string(msg1)+std::string(" ")+std::string(msg2))
  #else
    #define EXC_PRINT_2(msg1, msg2)                                                                                                \
    do {                                                                                                                        \
      std::cerr<<"file:"<<__FILE__<<", line:"<<__LINE__<<", function: "<<__PRETTY_FUNCTION__<<". "<<std::endl<<msg1<<" "<<msg2<<std::endl; \
      throw PrositAux::Exc(std::string(__PRETTY_FUNCTION__)+std::string(": ")+std::string(msg1)+std::string(" ")+std::string(msg2)); \
    } while(0)
  #endif

}

#endif
