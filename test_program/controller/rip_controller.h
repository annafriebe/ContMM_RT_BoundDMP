#ifndef RIP_CONTROLLER_H
#define RIP_CONTROLLER_H
#include <vector>
#include "rip_kf_estimator_sr_quanser.h"
#include "rot_pend_sensor_item.h"


class RIPController {
private:
	RIP_KFEstimatorSRQuanser& m_est;
	std::string m_est_filename;
	double m_voltage;
	ns::RotPendSensorItem m_si;
	double m_KP_pend;  
	double m_KD_pend;
	double m_KP_arm;
	double m_KD_arm;
	std::vector<ns::RotPendSensorItem> m_estimates;
	std::vector<int> m_missed_deadlines;
	int m_deadline_resolution;
public:
	RIPController(RIP_KFEstimatorSRQuanser& est, std::string est_filename,
		int deadline_resolution=500);

	void run(double sim_time_step, double calc_time_step, unsigned int n_steps,
		double rel_deadline);
	~RIPController();

};
#endif
