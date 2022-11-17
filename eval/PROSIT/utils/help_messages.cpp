/*!
 * @file    help_messages.cpp
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
#include "help_messages.hpp"

using namespace std;

namespace PrositAux {

  /// @brief Prints the help for the CLI probability solver.
  void help_cli() {

    cout << endl << "Usage:   ./cli_solver [options]" << endl
         << "Options:" << endl
         << "-c | --cyclic               Set the solving algorithm to: Cyclic Reduction." << endl
         << "-l | --logarithmic          Set the solving algorithm to: Latouche Reduction." << endl
         << "-a | --analytic             Set the solving algorithm to: Analytic Bound." << endl
         << "-o | --companion            Set the solving algorithm to: Companion Form." << endl
         << "-R | --reservation_based    Set the scheduling algorithm to: Resource Reservation." << endl
         << "-X | --fixed_priority       Set the scheduling algorithm to: Fixed Priority." << endl
         << "-T | --task_period          Period of the task." << endl
         << "-t | --server_period        Reservation period of the task." << endl
         << "-q | --budget               Budget of the task." << endl
         << "-C | --comp_time_path       Path to the file with the PMF of the computation times." << endl
         << "-Z | --inter_time_path      Path to the file with the PMF of the interarrival times." << endl
         << "-M | --trans_mat_path       Path to the file with the transition matrix of the system." << endl
         << "-d | --max_deadline         Maximum number of deadlines to analyze." << endl
         << "-g | --granularity          Granularity for resampling the probability distributions [1]." << endl
         << "-m | --num_modes            Number of modes in the system [1]." << endl
         << "-e | --epsilon              Tolerance for stopping the Latouche algorithm [1e-4]." << endl
         << "-i | --num_iterations       Maximum  number of  iterations  for the  Cyclic Reduction  or" << endl
         << "                            Latouche algorithm [100]." << endl
         << "-S | --shift_flag           Set the shift technique flag for Cyclic Reduction [false]." << endl
         << "-F | --qos_function         Type of the Quality of Service function [no QoS]." << endl
         << "-n | --qos_pmin             Lower bound of the QoS function [0.0]." << endl
         << "-x | --qos_pmax             Upper bound of the QoS function [0.0]." << endl
         << "-s | --qos_scale            Scale of the QoS function [0.0]." << endl
         << "-f | --qos_offset           Minimum value of the QoS function [0.0]." << endl
         << "-v | --verbose              Set the verbose mode to print online information [false]." << endl
         << "-h | --help                 Print this message." << endl << endl
         << "Example:" << endl 
         << "This command will  execute the cli-based Latouche  probability solver for a periodic task" << endl
         << "with  Markovian computation times  and scheduled using a  resource reservation scheduling" << endl
         << "algorithm.  The task period is:  100000 us, the reservation period is:  25000 us  and the" << endl
         << "budget is:  10000 us.  The PMFs of the computation times will be resampled by a factor of" << endl
         << "250 (ceil(sample/250)). The system is modeled as a 4-states Markov chain.                " << endl
         << "The probability distribution  of the computation  times for  each state is  stored in the" << endl
         << "files: tests/data/mctm/pmf_comp_time_state#.txt, while the transition matrix is stored in" << endl
         << "the file: tests/data/mctm/transition_matrix.txt. The program will present the probability" << endl
         << "of respecting  the deadline for the  first 12 multiples of the  reservation period.  Note" << endl
         << "that the last eleven parameters present a default value, hence they are optional." << endl << endl
         << "~$ ./cli_solver -R -l -T 100000 -t 25000 -q 10000 -g 250 -d 12 -m 4 -C tests/data/mctm/pmf_comp_time_state -M tests/data/mctm/transition_matrix.txt" 
         << endl
         << "./cli_solver -R -a -T 100000 -t 25000 -q 12500 -g 12500 -d 12 -C tests/data/iid/pmf_comp_times"
         << endl
         << "./cli_solver -R -o -T 100000 -t 25000 -q 12500 -g 250 -d 12 -C tests/data/iid/pmf_comp_times"
         << endl
         << "./cli_solver -R -c -T 100000 -t 25000 -q 12500 -g 250 -d 12 -C tests/data/iid/pmf_comp_times"
         << endl
         << "./cli_solver -X -d 4 -k 2 -P tests/data/fp/example3/parameters.txt -C tests/data/fp/example3/pmf_comp_time_task"
         << endl << endl;

    exit(0);

  }

  /// @brief Prints the help for the XML probability solver.
  void help_xml() {

    cout << endl << "Usage:   ./xml_solver [options]" << endl
         << "Options:" << endl
         << "-t | --taskset              Path to the xml file with the characteristics of the tasks." << endl
         << "-f | --spec_file            Path to the xml file with the parameters for the solver." << endl
         << "-v | --verbose              Set the verbose mode to print online information [false]." << endl
         << "-h | --help                 Print this message." << endl << endl
         << "Example:" << endl
         << "This command will  execute the xml-based probability solver  using the tasksets presented" << endl
         << "in the file:  examples/demo/taskset.xml, which are associated with the schedulers defined" << endl
         << "in the file: examples/demo/synthesis.xml. Note that the verbose option presents a default" << endl
         << "value, hence it is optional." << endl << endl
         << "~$ ./xml_solver -t examples/demo/taskset.xml -f examples/demo/synthesis.xml" 
         << endl << endl;

    exit(0);

  }

  /// @brief Prints the help for the autocorrelation program.
  void help_autocorr() {

    cout << endl << "Usage:   ./mctm_autocorrelation [options]" << endl
         << "Options:" << endl
         << "-o | --observation_path     Path to the file with the observation samples." << endl
         << "-l | --lags                 Number of time instants in the autocorrelation function [20]." << endl
         << "-d | --destination_path     Path  to the  file where  the autocorrelation  will be  saved" << endl
         << "                            [autocorrelation.txt]." << endl
         << "-h | --help                 Print this message." << endl << endl
         << "Example:" << endl 
         << "This command will execute the  autocorrelation function of the observations stored in the" << endl 
         << "file: tests/data/hmm/comp_times.txt. The function will be performed over 40 time instants" << endl
         << "and the results will be stored in the file: tests/data/hmm/autocorrelation.txt. Note that" << endl
         << "the last two parameters present a default value, hence they are optional." << endl << endl
         << "~$ ./mctm_autocorrelation -o tests/data/hmm/comp_times.txt -l 40 -d tests/data/hmm/autocorrelation/experimental.txt" 
         << endl << endl;

    exit(0);

  }

  /// @brief Prints the help for the HMM Learning problem.
  void help_learner() {

    cout << endl << "Usage:   ./mctm_learner [options]" << endl
         << "Options:" << endl
         << "-o | --observation_path     Path to the file with the observation samples." << endl
         << "-s | --num_states           Number of hidden states in the Markov chain." << endl
         << "-g | --granularity          Granularity for resampling the observations [1]." << endl
         << "-r | --num_trials           Number of random trials to choose the best solution [10]." << endl
         << "-t | --transition_path      Path  to the  file where  the transition  matrix will be saved" << endl
         << "                            [transition_matrix.txt]." << endl
         << "-e | --emission_path        Path  to the  files where  the emissions  matrix will be saved" << endl
         << "                            [pmf_state#.txt]." << endl
         << "-i | --num_iterations       Maximum  number of  iterations  for the  Baum-Welch  algorithm" << endl
         << "                            [1000]." << endl
         << "-p | --epsilon              Tolerance for stopping the Baum-Welch algorithm [1e-4]." << endl
         << "-h | --help                 Print this message." << endl << endl
         << "Example:" << endl 
         << "This command will execute the Baum-Welch algorithm for the HMM Learning problem. The input" << endl 
         << "observations are stored in the file:  tests/data/hmm/comp_times.txt. The algorithm defines" << endl
         << "a 5-states  Markov  chain  and the  observations will  be resampled  by a  factor  of 1000" << endl
         << "(ceil(sample/1000)).  The HMM parameters will be selected by the best likelihood among 150" << endl
         << "trials.  The destination path for the results is defined as:  tests/data/hmm/learner/. The" << endl
         << "stopping criteria for the Baum-Welch algorithm mantains  the default parameters (tolerance" << endl
         << "(error of 1e-4 or 1000 iterations).  Note that the last six  parameters present  a default" << endl
         << "value, hence they are optional." << endl
         << "This program  will present a table  with the likelihood  for each iteration along with the" << endl 
         << "best likelihood and the  corresponding best iteration.  Moreover, the program will  create" << endl 
         << "different text files containing the  probability transition matrix (transition_matrix.txt)" << endl 
         << "and the  probability distribution  for each state  (pmf_state#.txt).  These files  will be" << endl 
         << "created in the directory: tests/data/hmm/learner/" << endl << endl
         << "~$ ./mctm_learner -o tests/data/hmm/comp_times.txt -s 5 -r 150 -g 1000 -t tests/data/hmm/learner/transition_matrix.txt -e tests/data/hmm/learner/pmf_state"
         << endl << endl;

    exit(0);

  }

  /// @brief Prints the help for the HMM Decoding problem.
  void help_decoder() {

    cout << endl << "Usage:   ./mctm_decoder [options]" << endl
         << "Options:" << endl
         << "-o | --observation_path     Path to the file with the observation samples." << endl
         << "-t | --transition_path      Path to the file with the transition matrix." << endl
         << "-e | --emission_path        Path to the file with the emissions matrix." << endl
         << "-s | --num_states           Number of hidden states in the Markov chain." << endl
         << "-g | --granularity          Granularity for resampling the observations [1]." << endl
         << "-c | --states_path          Path to  the files  where the classified observations  will be" << endl
         << "                            saved [comp_time_state#.txt]." << endl
         << "-l | --significance         Significance level for the p-value [0.05]." << endl
         << "-h | --help                 Print this message." << endl << endl
         << "Example:" << endl 
         << "This command will execute  the Viterbi algorithm for the  HMM Decoding problem.  The input" << endl 
         << "observations are stored in the file:  tests/data/hmm/comp_times.txt; the transition matrix" << endl
         << "is stored in the file: tests/data/hmm/decoder/transition_matrix.txt, while the probability" << endl
         << "distribution for each state is stored in the files:  tests/data/hmm/decoder/pmf_state#.txt" << endl
         << "The algorithm defines a 5-states Markov chain and the observations will be resampled  by a" << endl
         << "factor of 1000  (ceil(sample/1000)).  The destination path  for the results is defined as:" << endl
         << "tests/data/hmm/decoder/.  The significance level for  the p-value is left with its default" << endl
         << "value (0.05).  Note that the last three parameters present a default value, hence they are" << endl
         << "optional." << endl
         << "This program  will present a table with the z-score and its corresponding p-value for each" << endl
         << "state in the Markov chain. Additionally, it shows whether the independence test was passed" << endl
         << "or not. Moreover, the program will create different text files containing the observations" << endl 
         << "classified by state. These files will be created in the directory: tests/data/hmm/decoder/" << endl << endl
         << "~$ ./mctm_decoder -o tests/data/hmm/comp_times.txt -s 5 -g 1000 -t tests/data/hmm/decoder/transition_matrix.txt -e tests/data/hmm/decoder/pmf_state -c tests/data/hmm/decoder/comp_time_state"
         << endl << endl;

    exit(0);

  }

  /// @brief Prints the help for the CBS simulator.
  void help_simul() {

    cout << endl << "Usage:   ./cbs_simulator [options]" << endl
         << "Options:" << endl
         << "-p | --is_pmf               Set the solving algorithm to: Cyclic Reduction." << endl
         << "-T | --task_period          Period of the task." << endl
         << "-t | --server_period        Reservation period of the task." << endl
         << "-q | --budget               Budget of the task." << endl
         << "-C | --comp_time_path       Path to the file with the PMF of the computation times." << endl
         << "-Z | --inter_time_path      Path to the file with the PMF of the interarrival times." << endl
         << "-M | --trans_mat_path       Path to the file with the transition matrix of the system." << endl
         << "-d | --max_deadline         Maximum number of deadlines to analyze." << endl
         << "-g | --granularity          Granularity for resampling the probability distributions [1]." << endl
         << "-m | --num_modes            Number of modes in the system [1]." << endl
         << "-l | --sim_length           Length of the simulation [100000]." << endl
         << "-c | --current_mode         Initial mode of the Markov chain [1]." << endl
         << "-r | --round_decimal        Number of decimals in the approximations [8]." << endl
         << "-e | --epsilon              Tolerance for the distributions [1e-4]." << endl
         << "-h | --help                 Print this message." << endl << endl
         << "Example:" << endl 
         << "This command will  execute the cli-based Latouche  probability solver for a periodic task" << endl
         << "with Markovian computation times. The task period is 100000 us, the reservation period is" << endl
         << "25000 us and the budget is 10000 us.  The PMFs of the computation times will be resampled" << endl
         << "by a factor of 250 (ceil(sample/250)).  The system is modeled as a 4-states Markov chain." << endl
         << "The probability distribution  of the computation  times for  each state is  stored in the" << endl
         << "files: tests/data/mctm/pmf_comp_time_state#.txt, while the transition matrix is stored in" << endl
         << "the file: tests/data/mctm/transition_matrix.txt. The program will present the probability" << endl
         << "of respecting  the deadline for the  first 12 multiples of the  reservation period.  Note" << endl
         << "that the last eleven parameters present a default value, hence they are optional." << endl << endl
         << "~$ ./cli_solver -l -T 100000 -t 25000 -q 10000 -g 250 -d 12 -m 4 -C tests/data/mctm/pmf_comp_time_state -M tests/data/mctm/transition_matrix.txt" 
         << endl << endl;

    exit(0);

  }

}
