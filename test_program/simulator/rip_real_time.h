#ifndef RIP_REAL_TIME_H
#define RIP_REAL_TIME_H

#include <memory>
#include <vector>
#include <mutex>
#include <string>
#include <atomic>
#include <random>
#include "quanser_rotary_inverted_pendulum_model.h"

namespace ns{
struct RotPendLogItem{
	double voltage;
	double pend_angle;
	double pend_ang_vel;
	double arm_angle;
	double arm_ang_vel;
	double end_time;
};


struct RotPendSensorItem{
	double pend_angle;
	double pend_ang_vel;
	int pend_angle_encoder;
	double arm_angle;
	double arm_ang_vel;
	int arm_angle_encoder;
	double time;	
};
}

class RIPRealTime {
private:
	QuanserRotaryInvertedPendulumModel& m_rip_model;
	std::vector<ns::RotPendLogItem> m_log_items;
	std::atomic<double> m_voltage;
	std::atomic<bool> m_startSignal;
	std::array<double, 4> m_state;
	double m_time;
	std::string m_log_file_name;
	std::unique_ptr<std::mt19937> m_gen_ptr;
	std::normal_distribution<> m_arm_sensor_noise;
	std::normal_distribution<> m_pend_sensor_noise;
	std::mutex m_lock;
	int m_log_resolution;

public:
	// constructor, starting from rest with arm and pendulum angle 0
	// TODO pass model with unique pointer, take ownership
	RIPRealTime(QuanserRotaryInvertedPendulumModel& rip_model, 
		std::string log_file_name, bool waitForController=true, int log_resolution = 100);
	// destructor, performing the logging to file
	~RIPRealTime();

	void runSimTime(double sim_time_step, double run_time_interval, unsigned int n_steps);
	ns::RotPendSensorItem getCurrentSensorItem();
	std::string getCurrentSensorItemJsonString();
	void startSimulation();
	void applyVoltage(double voltage);

};
#endif
