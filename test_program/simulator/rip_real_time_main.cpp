#include <pthread.h>
#include <sched.h>
#include <iostream>
#include "quanser_rotary_inverted_pendulum_model.h"
#include "rip_real_time.h"

using namespace std;

int main(int argc, char *argv[]) {
	int simfactor = 1;
	if (argc > 1){
		simfactor = stoi(argv[1]);
	}
 
	pthread_t this_thread = pthread_self();
  	struct sched_param params;
  	params.sched_priority = sched_get_priority_max(SCHED_FIFO);
  	int ret = pthread_setschedparam(this_thread, SCHED_FIFO, &params);
  	if (ret != 0) {
    	cerr << "Unable to set scheduling attributes." << endl;
  	}
  	else {
    	cerr << "Scheduling attributes properly set." << endl;
  	}
  	int policy = 0;
  	ret = pthread_getschedparam(this_thread, &policy, &params);
  	if (ret != 0) {
    	cerr << "Couldn't retrieve real-time scheduling paramers" << endl;
    	return 0;
  	}
  
  	// Check the correct policy was applied
  	if(policy != SCHED_FIFO) {
    	cerr << "Scheduling is NOT SCHED_FIFO!" << std::endl;
  	} else {
    	cerr << "SCHED_FIFO OK" << endl;
  	}
  	cpu_set_t set1;
  	CPU_ZERO(&set1);
  	CPU_SET(1, &set1);

  	ret = pthread_setaffinity_np(this_thread, sizeof(set1), &set1);
  	if (ret < 0){
    	cerr << "set affinity failed" << endl;
  	}

	QuanserRotaryInvertedPendulumModel rip(0.095 /*arm mass*/, 
	0.024 /* pend mass */, 0.085 /*arm len*/, 0.13 /*pend len*/, 
	0.0015 /* arm damping*/, 0.0005 /*pend damping*/, 
	0.042 /* motor torque constant*/, 0.042 /* back emf constant*/,
	8.4 /* motor torque resistance */, 1e-5, 3.14);
	RIPRealTime rip_rt(rip, "testRotaryRealTime.json");
	rip_rt.runSimTime(0.002, 0.002*simfactor, 50000);
}
