#ifndef ROT_PEND_SENSOR_ITEM_H
#define ROT_PEND_SENSOR_ITEM_H

#include "../nlohmann/json/single_include/nlohmann/json.hpp"

namespace ns{
struct RotPendSensorItem{
	double pend_angle;
	double pend_ang_vel;
	double arm_angle;
	double arm_ang_vel;
	double time;	
};

void to_json(nlohmann::json& j, const RotPendSensorItem& si);

}
#endif