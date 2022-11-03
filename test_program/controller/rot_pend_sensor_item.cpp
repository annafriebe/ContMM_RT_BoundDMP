#include "rot_pend_sensor_item.h"
#include "../nlohmann/json/single_include/nlohmann/json.hpp"

using json=nlohmann::json;
namespace ns{
void to_json(json& j, const RotPendSensorItem& si){
	j["pend angle"]= si.pend_angle;
	j["pend ang vel"] = si.pend_ang_vel;
	j["arm angle"] = si.arm_angle;
	j["arm ang vel"] = si.arm_ang_vel;
	j["time"] = si.time; 
}
	
}

