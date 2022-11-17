/*!
 * @file    cli_solver.cpp
 * 
 * @brief   A CLI solver for computation of stady state probabilities.
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
#include <getopt.h>

#include "../solver/companion_rr_probability_solver.hpp"
#include "../solver/analytic_rr_probability_solver.hpp"
#include "../solver/cr_qbd_rr_probability_solver.hpp"
#include "../solver/lt_qbd_rr_probability_solver.hpp"
#include "../solver/fp_probability_solver.hpp"
#include "../qos/quadratic_qos_function.hpp"
#include "../utils/auxiliary_functions.hpp"
#include "../tasks/rr_task_descriptor.hpp"
#include "../tasks/fp_task_descriptor.hpp"
#include "../qos/linear_qos_function.hpp"
#include "../tasks/fp_task_schedule.hpp"
#include "../utils/help_messages.hpp"

using namespace PrositCore;
using namespace PrositAux;
using namespace Eigen;
using namespace std;

static uint64_t T = 0;
static uint64_t Ts = 0;
static uint64_t Qs = 0;
static uint64_t granularity = 1;
static uint64_t max_deadline = 0;
static uint64_t num_modes = 1;
static uint64_t num_tasks = 1;
static int num_iter = 100;
static double eps = 1e-4;
static double qos_pmin = 0.0;
static double qos_pmax = 0.0;
static double qos_scale = 0.0;
static double qos_offset = 0.0;
static string algo = "cyclic";
static string path_comp_time;
static string path_inter_time;
static string path_trans_mat;
static string path_task_params;
static string path_results;
static string qos_type;
static bool verbose_flag = false;
static int cr_flag = 0;
static int logarithmic_flag = 0;
static int analytic_flag = 0;
static int companion_flag = 0;
static int shift_flag = 0;
static int fp_scheduler = 0;
static int rb_scheduler = 0;

/* Parse the command line arguments */
static int opts_parse(int argc, char *argv[]) {

  int opt;

  /* Short description for the arguments */
  static const char short_options[] = "T:t:q:g:d:m:e:i:C:Z:M:P:F:n:x:s:f:k:r:vclaoShRX";

  /* Long description for the arguments */
  static struct option long_options[] = {

    /* These options set a flag. */
    {"verbose",           no_argument, 0, 'v'},
    {"cyclic",            no_argument, 0, 'c'},
    {"logarithmic",       no_argument, 0, 'l'},
    {"analytic",          no_argument, 0, 'a'},
    {"companion",         no_argument, 0, 'o'},
    {"shift_flag",        no_argument, 0, 'S'},
    {"fixed_priority",    no_argument, 0, 'X'},
    {"reservation_based", no_argument, 0, 'R'},

    /* These options don't set a flag. We distinguish them by their indices */
    {"task_period",       required_argument, 0, 'T'},
    {"server_period",     required_argument, 0, 't'},
    {"budget",            required_argument, 0, 'q'},
    {"granularity",       required_argument, 0, 'g'},
    {"max_deadline",      required_argument, 0, 'd'},
    {"num_modes",         required_argument, 0, 'm'},
    {"num_tasks",         required_argument, 0, 'k'},
    {"epsilon",           required_argument, 0, 'e'},
    {"num_iterations",    required_argument, 0, 'i'},
    {"qos_function",      required_argument, 0, 'F'},
    {"qos_pmin",          required_argument, 0, 'n'},
    {"qos_pmax",          required_argument, 0, 'x'},
    {"qos_scale",         required_argument, 0, 's'},
    {"qos_offset",        required_argument, 0, 'f'},
    {"comp_time_path",    required_argument, 0, 'C'},
    {"inter_time_path",   required_argument, 0, 'Z'},
    {"trans_mat_path",    required_argument, 0, 'M'},
    {"task_params_path",  required_argument, 0, 'P'},
    {"results_path",      required_argument, 0, 'r'},
    {"help",              no_argument,       0, 'h'},
    {0, 0, 0, 0},

  };

  /* Iterate over the list of the arguments */
  while ((opt = getopt_long(argc, argv, short_options, long_options, 0)) != -1) {

    switch (opt) {

      /* Task period */
      case 'T':
        T = atoi(optarg);
        break;

      /* Reservation period */
      case 't':
        Ts = atoi(optarg);
        break;

      /* Budget */
      case 'q':
        Qs = atoi(optarg);
        break;

      /* Granularity */
      case 'g':
        granularity = atoi(optarg);
        break;

      /* Maximum deadline */
      case 'd':
        max_deadline = atoi(optarg);
        break;

      /* Number of MCTM modes */
      case 'm':
        num_modes = atoi(optarg);
        break;

      /* Epsilon error */
      case 'e':
        eps = atof(optarg);
        break;

      /* Maximum number of iterations */
      case 'i':
        num_iter = atoi(optarg);
        break;

      /* Path to the computation times PMFs */
      case 'C':
        path_comp_time = optarg;
        break;

      /* Path to the interarrival times PMF */
      case 'Z':
        path_inter_time = optarg;
        break;

      /* Path to the probability transition matrix */
      case 'M':
        path_trans_mat = optarg;
        break;

      /* Path to the parameters of the task */
      case 'P':
        path_task_params = optarg;
        break;

      /* Path to the computation times PMFs */
      case 'r':
        path_results = optarg;
        break;

      /* QoS function */
      case 'F':
        qos_type = optarg;
        break;

      /* Lower bound probability for the QoS function */
      case 'n':
        qos_pmin = atof(optarg);
        break;

      /* Upper bound probability for the QoS function */
      case 'x':
        qos_pmax = atof(optarg);
        break;

      /* Slope for the mapping of the QoS function */
      case 's':
        qos_scale = atof(optarg);
        break;

      /* Minimum value for the QoS function */
      case 'f':
        qos_offset = atof(optarg);
        break;

      /* Number of fixed-priority tasks */
      case 'k':
        num_tasks = atoi(optarg);
        break;

      /* Verbose flag */
      case 'v':
        verbose_flag = true;
        break;

      /* Cyclic Reduction flag */
      case 'c':
        cr_flag = 1;
        algo = "cyclic";
        break;

      /* Logarithmic Reduction flag */
      case 'l':
        logarithmic_flag = 1;
        algo = "logarithmic";
        break;

      /* Analitic flag */
      case 'a':
        analytic_flag = 1;
        algo = "analytic";
        break;

      /* Companion flag */
      case 'o':
        companion_flag = 1;
        algo = "companion";
        break;

      /* Shift flag */
      case 'S':
        shift_flag = 1;
        break;

      /* Set the resource reservation scheduler */
      case 'R':
        rb_scheduler = 1;
        break;

      /* Set the fixed priority scheduler */
      case 'X':
        fp_scheduler = 1;
        break;

      /* Print the help */
      case 'h':
        help_cli();
        break;

      default:
        EXC_PRINT("ERROR: Incorrect parameters during the parsing.");

    }

  }

  return optind;

}

/* Main program */
int main(int argc, char *argv[]) {

  /* Variables for the time analysis */
  long long t_start = 0, t_solution_start = 0, t_end = 0;

  /* Variables for the distributions ranges */
  uint64_t zmin = numeric_limits<uint64_t>::max();
  uint64_t zmax = numeric_limits<uint64_t>::min();

  try {

    /* Variable for the log files */
    ofstream results_file;

    /* Obtain the initial computation time */
    t_start = my_get_time();

    /* Parse the command line arguments */
    opts_parse(argc, argv);

    /* Check the selection of one scheduler */
    if ((fp_scheduler + rb_scheduler) != 1) {
      EXC_PRINT("ERROR: Ambiguous choice of the scheduler.");
    }

    /* The solver has to use the Reservation-Based Scheduling Algorithm */
    if (rb_scheduler) {

      /* Check the scheduling parameters */
      if ((Ts == 0) || (Qs == 0)) {
        EXC_PRINT("ERROR: The scheduling parameters have not been set.");
      }

      /* Set Cyclic Reduction as default solver */
      if ((analytic_flag + cr_flag + logarithmic_flag + companion_flag) == 0) {
        cr_flag = 1;
        algo = "cyclic";
      }

      /* Check the selection of one algorithm */
      if ((analytic_flag + cr_flag + logarithmic_flag + companion_flag) != 1) {
        EXC_PRINT("ERROR: Ambiguous choice of the algorithm.");
      }

      /* Check the maximum deadline to analyse */
      if (max_deadline <= 0) {
        EXC_PRINT("ERROR: The maximum deadline has not been properly set.");
      }
      else {
        
      }

      /* Check the correctness of the shift flag */
      if ((shift_flag != 0) && (cr_flag == 0) && (verbose_flag)) {
        cerr << "WARNING: The shift flag only makes sense for Cyclic Reduction." << endl;
      }

      /* Vector for the computation times */
      vector<unique_ptr<pmf>> comp_times;

      /* Iterate over the number of modes */
      for (uint32_t i = 1; i <= num_modes; i++) {
      //for (uint64_t i = 1; i <= num_modes; i++) {

        /* Variables for the distributions ranges */
        uint64_t cmin = numeric_limits<uint64_t>::max();
        uint64_t cmax = numeric_limits<uint64_t>::min();

        /* Variable for the name of the computation times */
        char comp_time_file[250];

        /* Build the string with the name */
        sprintf(comp_time_file, "%s%d.txt", path_comp_time.c_str(), i);
        //sprintf(comp_time_file, "%s%lu.txt", path_comp_time.c_str(), i);

        /* Check the PMF of the computation time */
        obtainRange(comp_time_file, cmin, cmax);
        //obtainRange(path_comp_time, cmin, cmax);

        /* Variable for the PMF of the computation time */
        unique_ptr<pmf> c(new pmf(cmin, cmax, eps));

        /* Load the PMF of the computation time */
        c->load(comp_time_file);
        //c->load(path_comp_time);

        /* Resample the computation times */
        //comp_times.push_back(move(c->resample(granularity, true)));
        comp_times.push_back(move(c));

      }

      /* The PMF for the interarrival time is passed as argument */
      if (!path_inter_time.empty()) {

        /* Check the PMF of the interarrival time */
        obtainRange(path_inter_time, zmin, zmax);

      }
      else {

        /* The task is periodic */
        if (T != 0) {

          /* Minimum and maximum indices equal to the task period */
          zmin = T;
          zmax = T;

        }
        else {
          EXC_PRINT("ERROR: The distribution of the interarrival time is missing.");
        }

      }

      /* Variable for the PMF of the interarrival time */
      unique_ptr<pmf> u(new pmf(zmin, zmax, eps));

      /* The PMF for the interarrival time is passed as argument */
      if (!path_inter_time.empty()) {

        /* Load the PMF of the interarrival time */
        u->load(path_inter_time);

      }
      else {

        /* Set the PMF of the interarrival time */
        u->set(0, 1.0);
        T = zmin;

      }

      /* Check whether the maximum deadline reaches the period */
      if (max_deadline < (T / Ts)) {
        max_deadline = T / Ts;
      }

      /* Variable for the transition matrix */
      MatrixXd transition_matrix = MatrixXd::Zero(num_modes, num_modes);
    
      /* Load the transition matrix */
      if (num_modes == 1) {
        transition_matrix(0, 0) = 1;
      }
      else {

        /* The probability transition matrix is passed as argument */
        if (!path_trans_mat.empty()) {
          loadMatrix(path_trans_mat, transition_matrix);
        }
        else {
          EXC_PRINT("ERROR: The probability transition matrix is missing.");
        }

      }

      /* Check the selection of the companion algorithm */
      if (companion_flag) {

        if (path_comp_time.empty()) {
          EXC_PRINT("ERROR: The distribution of the computation time is missing.");
        }

        /* TODO: Apply to aperiodic scenarios */
        if (T == 0) {
          EXC_PRINT("ERROR: The task period must be set for the Companion Form.");
        }

      }

      /* Obtain the final parsing time */
      t_solution_start = my_get_time();

      /* Create the task descriptor */
      ResourceReservationTaskDescriptor task_des("Tau1", num_modes, move(comp_times), move(u->resample(Ts, false)), transition_matrix, Qs, Ts, granularity, (max_deadline * Ts), algo);

      /* Set the parameters for the task descriptor */
      task_des.setDeadlineStep(Ts);
      task_des.setVerboseFlag(verbose_flag);

      /* Insert the deadline in the probability map */
      for (uint64_t i = 1; i <= max_deadline; i++) {
        task_des.insertDeadline(task_des.getDeadlineStep() * i);
      }

      /* Check whether the QoS function is present */
      if (!qos_type.empty()) {

        /* The QoS function is correct */
        if (strcmp(qos_type.c_str(), "linear") == 0) {

          /* Create the QoS function */
          unique_ptr<QoSFunction> qos(new LinearQoSFunction(qos_scale, qos_pmin, qos_pmax, qos_offset));
          task_des.q = move(qos);
          task_des.qos_type = qos_type;

        }
        else if (strcmp(qos_type.c_str(), "quadratic") == 0) {

          /* Create the QoS function */
          unique_ptr<QoSFunction> qos(new QuadraticQoSFunction(qos_scale, qos_pmin, qos_pmax));
          task_des.q = move(qos);
          task_des.qos_type = qos_type;

        }
        else {

          EXC_PRINT("ERROR: The QoS function type is not recognized.");

        }

      }

      /* Set the analytic solver */
      if (analytic_flag) {

        if ((num_modes == 1) && (task_des.isPeriodic())) {

          unique_ptr<ResourceReservationProbabilitySolver> ps(new AnalyticResourceReservationProbabilitySolver());
          task_des.setSolver(move(ps));
          task_des.computeProbability();

        }
        else {

          EXC_PRINT("ERROR: The analytic solver is only implemented for periodic single mode tasks.");

        }

      }
      else if (companion_flag) {

        if ((num_modes == 1) && (task_des.isPeriodic())) {
          unique_ptr<ResourceReservationProbabilitySolver> ps(new CompanionResourceReservationProbabilitySolver(/*max_deadline*/));
          task_des.setSolver(move(ps));
          task_des.computeProbability();
        }
        else {
          EXC_PRINT("ERROR: The companion solver is only implemented for periodic single mode tasks.");
        }

      }
      else if (cr_flag) {

        unique_ptr<CRQBDResourceReservationProbabilitySolver> ps(new CRQBDResourceReservationProbabilitySolver(shift_flag ? true : false, num_iter));
        task_des.setSolver(move(ps));
        task_des.computeProbability();

      }
      else {
      
        unique_ptr<LatoucheQBDResourceReservationProbabilitySolver> ps(new LatoucheQBDResourceReservationProbabilitySolver(eps, num_iter));
        task_des.setSolver(move(ps));
        task_des.computeProbability();

      }

      /* Present the header for the results */
      cout << "===========================================================================" << endl;
      cout << "=                                 Results                                 =" << endl;
      cout << "===========================================================================" << endl;
      cout << setw(17) << "Deadline" << setw(26) << "Probability" << setw(23) << "Quality" << endl;

      /* The name of the result file is passed as argument */
      if (!path_results.empty()) {

        /* Open the log file */
        results_file.open(path_results);

        /* Set the initial value to zero */
        results_file << "0 0.00000000000000000000" << endl;

      }

      /* Present the probability of respecting the deadlines */
      for (uint64_t i = 1; i <= max_deadline; i++) {

        double quality;
        double probability = task_des.getProbability(task_des.getDeadlineStep() * i);

        /* Check whether the QoS function exists */
        if (!(task_des.q)) {

          /* No QoS function */
          quality = 0.0;

        }
        else {

          /* QoS function */
          quality = task_des.q->evaluateQuality(probability);

        }

        cout << setw(17) << fixed << task_des.getDeadlineStep() * i 
             << setw(31) << setprecision(20) << probability 
             << setw(18) << setprecision(4)  << quality
             << setprecision(2) << endl;

        /* The name of the result file is passed as argument */
        if (!path_results.empty()) {

          results_file << fixed << task_des.getDeadlineStep() * i << " "
                       << setprecision(20) << probability << endl;

        }

      }

      /* The name of the result file is passed as argument */
      if (!path_results.empty()) {

        /* Close the log file */
        results_file.close();

      }

      /* Obtain the final computation time */
      t_end = my_get_time();

      /* Present the final results */
      cout << "===========================================================================" << endl;
      cout << setw(9)  << "Name:"     << setw(26) << task_des.getName() << endl;
      cout << setw(11) << "Budget:" << setw(24) << task_des.getBudget() << endl;
      cout << setw(14) << "Bandwidth:" << setw(21) << double(task_des.getBudget()) / double(task_des.getServerPeriod()) << endl;
      cout << "===========================================================================" << endl;
      cout << "=                            Computation  Time                            =" << endl;
      cout << "===========================================================================" << endl;
      cout << setw(17) << "Parsing time:"  << setw(18) << t_solution_start - t_start << " us" << endl;
      cout << setw(18) << "Solution time:" << setw(17) << t_end - t_solution_start   << " us" << endl;
      cout << setw(15) << "Total time:"    << setw(20) << t_end - t_start            << " us" << endl;
      cout << "===========================================================================" << endl;

      /* Delete the object */
      u.reset();

    }
    else if (fp_scheduler) {

      FixedPriorityTaskSet taskset;

      /* Variable for the task parameters */
      MatrixXd task_parameters = MatrixXd::Zero(num_tasks, 3);

      /* Load the task parameters */
      if (!path_task_params.empty()) {
        loadMatrix(path_task_params, task_parameters);
      }
      else {
        EXC_PRINT("ERROR: The task parameters are missing.");
      }

      /* Iterate over the task set */
      for (uint32_t i = 0; i < num_tasks; i++) {
      //for (uint64_t i = 0; i < num_tasks; i++) {

        /* Obtain the period of the task */
        T = task_parameters(i, 0);

        /* Variable for the PMF of the interarrival time */
        unique_ptr<pmf> u(new pmf(T, T, eps));

        /* Set the PMF of the interarrival time */
        u->set(0, 1.0);

        /* Variables for the distributions ranges */
        uint64_t cmin = numeric_limits<uint64_t>::max();
        uint64_t cmax = numeric_limits<uint64_t>::min();

        /* Variable for the name of the computation times */
        char comp_time_file[250];

        /* Build the string with the name */
        sprintf(comp_time_file, "%s%d.txt", path_comp_time.c_str(), i + 1);
        //sprintf(comp_time_file, "%s%lu.txt", path_comp_time.c_str(), i + 1);

        /* Check the PMF of the computation time */
        obtainRange(comp_time_file, cmin, cmax);

        /* Variable for the PMF of the computation time */
        unique_ptr<pmf> c(new pmf(cmin, cmax, eps));

        /* Load the PMF of the computation time */
        c->load(comp_time_file);

        /* Variable for the name of the task */
        char task_name[20];
        //char task_name[6];

        /* Build the string with the name */
        sprintf(task_name, "Tau%d", i + 1);
        //sprintf(task_name, "Tau%lu", i + 1);

        /* Create the task descriptor */
        FixedPriorityTaskDescriptor *task_des = new FixedPriorityTaskDescriptor(task_name, move(c->resample(granularity, true)), move(u), 20, task_parameters(i, 1), task_parameters(i, 2));

        /* Insert the deadline in the probability map */
        for (uint64_t i = 1; i <= max_deadline; i++) {
          task_des->insertDeadline(task_des->getDeadlineStep() * i);
        }

        /* Check whether the QoS function is present */
        if (!qos_type.empty()) {

          /* The QoS function is correct */
          if (strcmp(qos_type.c_str(), "linear") == 0) {

            /* Create the QoS function */
            unique_ptr<QoSFunction> qos(new LinearQoSFunction(qos_scale, qos_pmin, qos_pmax, qos_offset));
            task_des->q = move(qos);
            task_des->qos_type = qos_type;

          }
          else if (strcmp(qos_type.c_str(), "quadratic") == 0) {

            /* Create the QoS function */
            unique_ptr<QoSFunction> qos(new QuadraticQoSFunction(qos_scale, qos_pmin, qos_pmax));
            task_des->q = move(qos);
            task_des->qos_type = qos_type;

          }
          else {

            EXC_PRINT("ERROR: The QoS function type is not recognized.");

          }

        }

        /* Include the task in the task set */
        taskset.push_back(move(task_des));

      }

      /* Obtain the final parsing time */
      t_solution_start = my_get_time();

      unique_ptr<FixedPriorityTaskSchedule> task_sched(new FixedPriorityTaskSchedule("Schedule1", move(taskset)));

      /* Set the parameters for the task schedule */
      task_sched->setVerboseFlag(verbose_flag);

      unique_ptr<FixedPriorityProbabilitySolver> ps(new FixedPriorityProbabilitySolver(eps, num_iter));
      task_sched->setSolver(move(ps));
      task_sched->computeProbability();

      /* Present the header for the results */
      cout << "===========================================================================" << endl;
      cout << "=                                 Results                                 =" << endl;
      cout << "===========================================================================" << endl;
      cout << setw(9)  << "Name:"     << setw(26) << task_sched->getName() << endl;
      cout << "===========================================================================" << endl;

      /* Iterate over the task */
      for (size_t i = 0; i < num_tasks; i++) {
      
        cout << "Task" << setw(13) << "Deadline" << setw(26) << "Probability" << setw(23) << "Quality" << endl;

        /* Present the probability of respecting the deadlines */
        for (uint64_t j = 1; j <= max_deadline; j++) {

          double quality;// = 0.0;
          double probability = task_sched->getProbability(((task_sched->getTaskSet())[i])->getDeadlineStep() * j, i);

          /* Check whether the QoS function exists */
          if (!(((task_sched->getTaskSet())[i])->q)) {

            /* No QoS function */
            quality = 0.0;

          }
          else {

            /* QoS function */
            quality = (((task_sched->getTaskSet())[i])->q)->evaluateQuality(probability);

          }

          cout << ((task_sched->getTaskSet())[i])->getName() 
               << setw(13) << fixed << ((task_sched->getTaskSet())[i])->getDeadlineStep() * j 
               << setw(31) << setprecision(20) << probability 
               << setw(18) << setprecision(4)  << quality
               << setprecision(2) << endl;

        }

        cout << "===========================================================================" << endl;

      }

      /* Obtain the final computation time */
      t_end = my_get_time();

      /* Present the final results */
      cout << "=                            Computation  Time                            =" << endl;
      cout << "===========================================================================" << endl;
      cout << setw(17) << "Parsing time:"  << setw(18) << t_solution_start - t_start << " us" << endl;
      cout << setw(18) << "Solution time:" << setw(17) << t_end - t_solution_start   << " us" << endl;
      cout << setw(15) << "Total time:"    << setw(20) << t_end - t_start            << " us" << endl;
      cout << "===========================================================================" << endl;

    }
    else {

      EXC_PRINT("ERROR: The scheduling algorithm is not recognized.");

    }

  }
  catch (Exc &e) {

    cerr << "Exception caught" << endl;
    e.what();

  }

  return 0;

}
