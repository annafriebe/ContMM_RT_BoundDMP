#include <fstream>
#include <iostream>
#include <chrono>
#include <thread>
#include <numbers>
#include <boost/bind/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/asio.hpp>
#include "../nlohmann/json/single_include/nlohmann/json.hpp"
#include "rip_controller.h"

using namespace std;
using json=nlohmann::json;
using boost::asio::ip::tcp;



namespace{
const unsigned short portn = 5020;
string port="5020";	

double angleFromEncoder(int encoderValue){
	// angle in radians
	return encoderValue * std::numbers::pi / 1024;
}

}


// TODO client side


RIPController::RIPController(RIP_KFEstimatorSRQuanser& est, std::string est_filename,
	int deadline_resolution): 
m_est(est), 
m_est_filename(est_filename),
m_voltage(0), 
m_KP_pend(30), 
m_KD_pend(2.5),
m_KP_arm(2),
m_KD_arm(2),
m_deadline_resolution(deadline_resolution){
	m_si.pend_angle = std::numbers::pi;
	m_si.pend_ang_vel = 0;
	m_si.arm_angle = 0;
	m_si.arm_ang_vel = 0;
}


void RIPController::run(double sim_time_step, double calc_time_step, unsigned int n_steps, double rel_deadline) {
	boost::asio::io_context io_context;

  	tcp::resolver resolver(io_context);

  	tcp::resolver::query query("localhost", port);
  	tcp::resolver::results_type endpoints =
    	resolver.resolve(query);
	m_missed_deadlines.resize(n_steps/m_deadline_resolution+1, 0);


  	auto start = chrono::steady_clock::now();
	
	auto nanosecond_duration = chrono::nanoseconds((int(calc_time_step * 1e9)));
  	auto deadline = start + rel_deadline*nanosecond_duration;
	tcp::socket socket(io_context);
	boost::asio::connect(socket, endpoints);

	for (unsigned int i=0; i<n_steps; ++i){
		// connect to the server and get sensor data
		// assuming the data fits in this buffer
		array<char, 512> buf;
		boost::system::error_code error;

		size_t len = socket.read_some(boost::asio::buffer(buf), error);
   		string sensor_string(&buf[0], len);
   		json j = json::parse(sensor_string);
		m_si.pend_angle = angleFromEncoder(j["pend angle enc"]);
		m_si.arm_angle = angleFromEncoder(j["arm angle enc"]);
   		Eigen::Vector2d y(m_si.arm_angle, m_si.pend_angle);
 		m_est.update(m_voltage, y, sim_time_step);
		double pend_angle_est = m_est.getPendAngleEstimate();
		double pend_angvel_est = m_est.getPendAngVelEstimate();
		double arm_angle_est = m_est.getArmAngleEstimate();
		double arm_angvel_est = m_est.getArmAngVelEstimate();
		m_voltage = m_KP_pend*sin(pend_angle_est) - m_KD_pend*pend_angvel_est
			+ m_KP_arm*sin(arm_angle_est) + m_KD_arm * arm_angvel_est;

		string voltage_str = to_string(m_voltage);
		boost::asio::write(socket, boost::asio::buffer(voltage_str));

		if (chrono::steady_clock::now() > deadline){
			m_missed_deadlines[i/m_deadline_resolution] +=1;
		}
		this_thread::sleep_until(start + nanosecond_duration);
		start += nanosecond_duration;
		deadline += nanosecond_duration;
	}
}


RIPController::~RIPController(){
	json j;
	json est_items = m_estimates;
	j["estimates"] = est_items;
	std::map<unsigned int, int> missed_deadline_ind;
	int sum_missed_deadlines = 0;
	for (unsigned int i = 0; i<m_missed_deadlines.size(); ++i){
		if (m_missed_deadlines[i] > 0){
			missed_deadline_ind  [i*m_deadline_resolution] = m_missed_deadlines[i];
			sum_missed_deadlines += m_missed_deadlines[i];
		}
	}
	json deadlines(missed_deadline_ind);
	j["missed deadlines"] = deadlines; 
	j["sum missed deadlines"] = sum_missed_deadlines;
	std::ofstream f(m_est_filename);
	f << setw(4) << j << std::endl;
}

