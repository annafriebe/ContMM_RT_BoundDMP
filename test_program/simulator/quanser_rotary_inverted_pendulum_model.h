#ifndef QUANSER_ROT_INV_PEND_MODEL_H
#define QUANSER_ROT_INV_PEND_MODEL_H

#include <array>

// Using Räsänen and Pyrhönen, 2019
class QuanserRotaryInvertedPendulumModel {
private:
	// SI units
	double m_pend_mass;
	double m_arm_length;
	double m_pend_length;
	double m_arm_damping;
	double m_pend_damping;
	double m_arm_inertia_moment;
	double m_pend_inertia_moment;
	double m_motor_torque_constant;
	double m_motor_back_EMF_constant;
	double m_motor_armature_resistance;
	double m_time;
	double m_voltage;
	// approximate/ starting time delta
	double m_time_delta;
	// state vector is arm angle, pendulum angle, arm anglular velocity, pendulum angular velocity
	std::array<double, 4> m_state;

public:
	// constructor, starting from rest with arm angle 0 and pendulum angle 0
	QuanserRotaryInvertedPendulumModel(double arm_mass, double pend_mass, 
		double arm_length, double pend_length, double arm_damping, 
		double pend_damping, double motor_torque_constant, 
		double motor_back_EMF_constant, double motor_armature_resistance, 
		double time_delta, double init_pend_angle);

	// operator for call to integrator of the differential equation
	void operator() (const std::array<double, 4>& state, std::array<double, 4>& state_derivative, const double /*time*/ );

	void simulateWithVoltage(double voltage, double sim_time_duration, std::array<double, 4>& result_state, double& stop_time);

	double getHorizontalArmLength() const;
	double getPendulumArmLength() const;


};
#endif
