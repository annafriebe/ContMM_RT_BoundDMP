#ifndef RIP_KF_ESTIMATOR_SR_QUANSER_H
#define RIP_KF_ESTIMATOR_SR_QUANSER_H

#include <Eigen/Dense>


class RIP_KFEstimatorSRQuanser {
private:

	// state vector is arm angle, pendulum angle, arm anglular velocity, pendulum angular velocity
	Eigen::MatrixXd m_state;
	Eigen::MatrixXd m_H; 
	Eigen::Matrix2d m_R_SR; // measurement error cov
	Eigen::Matrix4d m_P_SR;
	Eigen::Matrix4d m_F;
	Eigen::Matrix4d m_Fderiv; 
	Eigen::Matrix4d m_Q_SR; // state error cov square root
	Eigen::MatrixXd m_G; // control signal to state
	Eigen::MatrixXd m_Gderiv;


public:
	RIP_KFEstimatorSRQuanser();

	void update(double voltage, Eigen::Vector2d y, double time_step);

	double getPendAngleEstimate() const;
	double getPendAngVelEstimate() const;
	double getArmAngleEstimate() const;
	double getArmAngVelEstimate() const;

};
#endif