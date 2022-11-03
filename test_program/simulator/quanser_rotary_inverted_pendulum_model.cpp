#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <cstdio>
#include "quanser_rotary_inverted_pendulum_model.h"


using namespace std;
using namespace boost::numeric::odeint;

namespace {
	const double earth_gravity_const = 9.81;
}

QuanserRotaryInvertedPendulumModel::QuanserRotaryInvertedPendulumModel(
	double arm_mass, double pend_mass, double arm_length, double pend_length, 
	double arm_damping, double pend_damping, double motor_torque_constant, 
	double motor_back_EMF_constant, double motor_armature_resistance, 
	double time_delta, double init_pend_angle): 
m_pend_mass(pend_mass),
m_arm_length(arm_length), 
m_pend_length(pend_length),
m_arm_damping(arm_damping),
m_pend_damping(pend_damping),
m_motor_torque_constant(motor_torque_constant),
m_motor_back_EMF_constant(motor_back_EMF_constant),
m_motor_armature_resistance(motor_armature_resistance),
m_time(0),
m_voltage(0),
m_time_delta(time_delta),
m_state({0, init_pend_angle, 0, 0}){
	m_arm_inertia_moment = arm_mass * pow(arm_length, 2)/12;
	m_pend_inertia_moment = pend_mass * pow(pend_length, 2)/12;
}



void QuanserRotaryInvertedPendulumModel::operator()
(const array<double, 4>& state, array<double, 4>& state_derivative, const double /*time*/ ){
	state_derivative[0] = state[2];
	state_derivative[1] = state[3];
	double denom1 = pow(0.5 * m_pend_mass * m_pend_length * m_arm_length * 
		cos(state[1]), 2) + 
	(m_pend_mass * pow(m_arm_length, 2) + m_arm_inertia_moment + 
		m_pend_mass*pow(0.5*m_pend_length*sin(state[1]), 2))*
	(m_pend_inertia_moment + 0.5*m_pend_mass*pow(m_pend_length, 2));
	double sd2term1 = -(m_pend_inertia_moment + 
		0.5*m_pend_mass*pow(m_pend_length, 2))*0.5*m_pend_mass*
	pow(m_pend_length, 2)*sin(state[1])*cos(state[1])*state[2]*state[3];
	double sd2term2 = - (m_pend_inertia_moment + 
		0.5*m_pend_mass*pow(m_pend_length, 2))*0.5*m_pend_mass*m_pend_length*
	m_arm_length*sin(state[1])*pow(state[3],2);
	double tau = m_motor_back_EMF_constant* (m_voltage - 
		m_motor_torque_constant*state[2])/m_motor_armature_resistance;
	double sd2term3 = (m_pend_inertia_moment + 0.5*m_pend_mass*
		pow(m_pend_length, 2))*(tau - m_arm_damping*state[2]);
	double sd2term4 = (0.5*m_pend_mass*m_pend_length*m_arm_length*cos(state[1]))*
	(0.25*m_pend_mass*pow(m_pend_length, 2)*cos(state[1])*sin(state[1])-
		0.5*m_pend_mass*m_pend_length*earth_gravity_const*sin(state[1]) -
		m_pend_damping*state[3]);
	state_derivative[2] = (sd2term1 + sd2term2 + sd2term3 + sd2term4)/denom1;
	double denom2 = m_pend_inertia_moment + m_pend_mass*
	pow(0.5*m_pend_length, 2);
	double sd3term1 = m_pend_mass*pow(0.5*m_pend_length, 2)*
	cos(state[1])*sin(state[1])*pow(state[3],2);
	double sd3term2 = -0.5*m_pend_mass*m_pend_length*m_arm_length*
	cos(state[1])*state_derivative[2];
	double sd3term3 = -0.5*m_pend_mass*m_pend_length*earth_gravity_const*
	sin(state[1]);
	double sd3term4 = -m_pend_damping*state[3];
	state_derivative[3] = (sd3term1 + sd3term2 + sd3term3 + sd3term4)/denom2;
}

void QuanserRotaryInvertedPendulumModel::simulateWithVoltage(double voltage, double sim_time_duration, std::array<double, 4>& result_state, double& stop_time) {
	m_voltage = voltage;
	integrate(*this, m_state, m_time, m_time + sim_time_duration, m_time_delta);
	m_time += sim_time_duration;
	for (unsigned int i=0; i<4; ++i){
		result_state[i] = m_state[i];
	}
	stop_time = m_time;
}

double QuanserRotaryInvertedPendulumModel::getHorizontalArmLength() const{
	return m_arm_length;
}

double QuanserRotaryInvertedPendulumModel::getPendulumArmLength() const{
	return m_pend_length;
}

