#include "rip_kf_estimator_sr_quanser.h"
#include <iostream>
#include <numbers>

using namespace Eigen;

namespace {
	const double earth_gravity_const = 9.81;
}



RIP_KFEstimatorSRQuanser::RIP_KFEstimatorSRQuanser()
{
	// initial estimate 
	m_state = MatrixXd::Zero(4,1);
	m_H = MatrixXd(2,4);
	m_H << MatrixXd::Identity(2,2), MatrixXd::Zero(2,2);

	m_R_SR << MatrixXd::Identity(2,2);
	m_R_SR(0,0) *= sqrt(0.02);
	m_R_SR(1,1) *= sqrt(0.02);

	m_P_SR << MatrixXd::Identity(4,4);
	m_P_SR(0,0) *= sqrt(0.02);
	m_P_SR(1,1) *= sqrt(0.02);
	m_P_SR(2,2) *= sqrt(0.2);
	m_P_SR(3,3) *= sqrt(0.2);

	Matrix4d Q;
	Q << MatrixXd::Identity(4,4);
	Q(0,0) *= 0.02;
	Q(1,1) *= 0.02;
	Q(2,2) *= 0.5;
	Q(3,3) *= 0.5;
	m_Q_SR = Q.llt().matrixL();

	m_Fderiv = MatrixXd::Zero(4,4);
	m_Fderiv(0,2) = 1;
	m_Fderiv(1,3) = 1;
	m_Fderiv(2,1) = -41.6;
	m_Fderiv(2,2) = -4.16;
	m_Fderiv(2,3) = 1.37;
	m_Fderiv(3,1) = 72.4;
	m_Fderiv(3,2) = -4.11;
	m_Fderiv(3,3) = -2.4;

	m_Gderiv = MatrixXd::Zero(4,1);
	m_Gderiv(2,0) = 13.9;
	m_Gderiv(3,0) = 13.7;

}


void RIP_KFEstimatorSRQuanser::update(double voltage, Eigen::Vector2d y, double time_step) {

	m_F = MatrixXd::Identity(4,4);
	m_F += time_step*m_Fderiv;

	m_G = time_step*m_Gderiv;

	MatrixXd P = m_Q_SR*m_Q_SR.transpose() + m_F*m_P_SR*m_P_SR.transpose()*m_F.transpose();
	// time update
	m_state = m_F*m_state + m_G * voltage;

	Matrix2d S = m_R_SR*m_R_SR + m_H * P * m_H.transpose();
	MatrixXd K = P * m_H.transpose() * S.inverse();

	// measurement update
	MatrixXd yMatrix = MatrixXd::Zero(2,1);
	yMatrix(0,0) = y(0);
	yMatrix(1,0) = y(1) - std::numbers::pi;

	Vector2d innovation = yMatrix - m_H*m_state;
	m_state = m_state + K*innovation;

	P = (MatrixXd::Identity(4,4) - K*m_H)*P;
	m_P_SR = P.llt().matrixL();

}

double RIP_KFEstimatorSRQuanser::getPendAngleEstimate() const {
	return m_state(1) + std::numbers::pi;
}
                                                                            
double RIP_KFEstimatorSRQuanser::getPendAngVelEstimate() const {
	return m_state(3);
}

double RIP_KFEstimatorSRQuanser::getArmAngleEstimate() const{
	return m_state(0);	
}

double RIP_KFEstimatorSRQuanser::getArmAngVelEstimate() const{
	return m_state(2);	
}
