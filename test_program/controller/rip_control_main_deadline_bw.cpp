#include <pthread.h>
#include <sched.h>
#include <iostream>
#include "rip_controller.h"
#include "rip_kf_estimator_sr_quanser.h"
#include "dl_syscalls.h"

using namespace std;

int main(int argc, char *argv[]) {
    int runtime = 70000;
    int period = 500000;
    int n = 4;
    int k=8;
    string filename = "estRotaryKalman.json";
    int simfactor = 1;
    if (argc < 7){
      cerr << "run with arguments: runtime, server period, n, k, output filename, simfactor" << endl;
      return -1;
    }
    runtime = stoi(argv[1]);
    period = stoi(argv[2]);
    n = stoi(argv[3]);
    k = stoi(argv[4]);
    filename = argv[5];
    simfactor = stoi(argv[6]);

  	pthread_t this_thread = pthread_self();
  	cpu_set_t set2;
  	CPU_ZERO(&set2);
  	CPU_SET(2, &set2);

  	struct sched_attr attr;
    memset(&attr, 0, sizeof(attr));
  	attr.size = sizeof(attr);
  	attr.sched_policy = SCHED_DEADLINE;
  	// Test high bandwidth
  	attr.sched_runtime =   runtime;
  	attr.sched_period =   period;
  	attr.sched_deadline = period;
  	attr.sched_flags = 0;
    int res = sched_setattr(0, &attr, 0);
    if (res < 0){
      int errsv = errno;
      cerr << "sched_setattr failed: " << endl;
      if (errsv == EINVAL){
        cerr << "EINVAL" << endl;
      }
      if (errsv == ESRCH){
        cerr << "ESRCH" << endl;
      }
      if (errsv == E2BIG){
        cerr << "E2BIG" << endl;
      }
      if (errsv == EBUSY){
        cerr << "EBUSY" << endl;
      }
      if (errsv == EPERM){
        cerr << "EPERM" << endl;
      }
    }
  res = sched_getattr(0, &attr, sizeof(attr), 0);
  if (res < 0){
    cerr << "sched_getattr failed" << endl;
  }
  if (attr.sched_policy != SCHED_DEADLINE){
    cerr << "SCHED_DEADLINE not set" << endl;
  } 
  res = pthread_setaffinity_np(this_thread, sizeof(set2), &set2);
  if (res < 0){
    cerr << "set affinity failed" << endl;
  }


 	RIP_KFEstimatorSRQuanser test;

  double rel_deadline = k/n;
	RIPController rip_controller(test, filename);
	rip_controller.run(0.002, 0.002*simfactor, 50000, rel_deadline);
}
