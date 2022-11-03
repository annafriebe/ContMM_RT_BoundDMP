#include <fstream>
#include <iostream>
#include <chrono>
#include <thread>
#include <numbers>
#include <boost/bind/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/asio.hpp>
#include <boost/math/special_functions/round.hpp>
#include "../nlohmann/json/single_include/nlohmann/json.hpp"
#include "rip_real_time.h"

using namespace std;
using json=nlohmann::json;
using boost::asio::ip::tcp;


namespace ns{
void to_json(json& j, const RotPendLogItem& li){
	j["voltage"] = li.voltage;
	j["pend angle"]= li.pend_angle;
	j["pend ang vel"] = li.pend_ang_vel;
	j["arm angle"] = li.arm_angle;
	j["arm ang vel"] = li.arm_ang_vel;
	j["end time"] = li.end_time; 
}

void to_json(json& j, const RotPendSensorItem& si){
	j["pend angle"]= si.pend_angle;
	j["pend ang vel"] = si.pend_ang_vel;
	j["pend angle enc"] = si.pend_angle_encoder;
	j["arm angle"] = si.arm_angle;
	j["arm ang vel"] = si.arm_ang_vel;
	j["arm angle enc"] = si.arm_angle_encoder;
	j["time"] = si.time; 
}

}

namespace{
	const unsigned short portn = 5020;

	int encoderFromAngle(double angleRadians){
		int encoder = boost::math::iround(
			1024 * angleRadians/std::numbers::pi);
		return encoder;
	}

	class tcp_connection
  		: public boost::enable_shared_from_this<tcp_connection>
	{
		public:
  		typedef boost::shared_ptr<tcp_connection> pointer;

  		static pointer create(boost::asio::io_context& io_context, RIPRealTime& rip_rt)
  		{
    		return pointer(new tcp_connection(io_context, rip_rt));
  		}

  		tcp::socket& socket()
  		{
    		return socket_;
  		}

  		void start()
  		{
  			/* Get current sensor json string */
  			message_ = rip_rt_.getCurrentSensorItemJsonString();

			boost::asio::async_write(socket_, boost::asio::buffer(message_),
				boost::bind(&tcp_connection::handle_write, shared_from_this(),
				boost::asio::placeholders::error,
				boost::asio::placeholders::bytes_transferred));
  		}

		private:
 		tcp_connection(boost::asio::io_context& io_context, RIPRealTime& rip_rt)
    		: socket_(io_context), rip_rt_(rip_rt)
  		{
  		}

  		void handle_read(const boost::system::error_code& error,
			size_t bytes_transferred)
  		{
			if (!error)
			{
				string voltage_string(&voltage_buf_[0], bytes_transferred);
				rip_rt_.applyVoltage(stod(voltage_string));
				/* Get current sensor json string */
	  			message_ = rip_rt_.getCurrentSensorItemJsonString();

				boost::asio::async_write(socket_,
					boost::asio::buffer(message_),
					boost::bind(&tcp_connection::handle_write, shared_from_this(),
				boost::asio::placeholders::error,
					boost::asio::placeholders::bytes_transferred));
			}
			else
			{
				// Ignoring error
			}
  		}


  		void handle_write(const boost::system::error_code& error,
      		size_t /*bytes_transferred*/)
  		{
			/* After write sensor data expect voltage command */
			if (!error)
			{
				socket_.async_read_some(boost::asio::buffer(voltage_buf_),
					boost::bind(&tcp_connection::handle_read, shared_from_this(),
				boost::asio::placeholders::error,
				boost::asio::placeholders::bytes_transferred));
				rip_rt_.startSimulation();

			}
			else
			{
				// Ignoring error
			}
  		}

  		tcp::socket socket_;
  		std::string message_;
  		RIPRealTime& rip_rt_;
  		array<char, 256> voltage_buf_;

	};
	class tcp_server
	{
		public:
  		tcp_server(boost::asio::io_context& io_context, RIPRealTime& rip_rt)
		: io_context_(io_context),
		acceptor_(io_context, tcp::endpoint(tcp::v4(), portn)), rip_rt_(rip_rt)
  		{
			start_accept();
  		}

		private:
  		void start_accept()
  		{
			tcp_connection::pointer new_connection =
				tcp_connection::create(io_context_, rip_rt_);

			acceptor_.async_accept(new_connection->socket(),
				boost::bind(&tcp_server::handle_accept, this, new_connection,
				boost::asio::placeholders::error));
  		}

  		void handle_accept(tcp_connection::pointer new_connection,
      		const boost::system::error_code& error)
  		{
			if (!error)
			{
				new_connection->start();
			}
			start_accept();
  		}

  		boost::asio::io_context& io_context_;
  		tcp::acceptor acceptor_;
  		RIPRealTime& rip_rt_;
	};

	void run_server_thread(RIPRealTime& rip_rt, boost::asio::io_context& io_context){
		try
  		{
			tcp_server server(io_context, rip_rt);
			io_context.run();
  		}
  		catch (std::exception& e)
  		{
			std::cerr << e.what() << std::endl;
  		}
	}

}


RIPRealTime::RIPRealTime(
	QuanserRotaryInvertedPendulumModel& rip_model,
	std::string log_file_name, bool waitForController, int log_resolution): 
	m_rip_model(rip_model),
	m_log_file_name(log_file_name),
	m_arm_sensor_noise(0,0.001),
	m_pend_sensor_noise(0,0.001),
	m_log_resolution(log_resolution)
	{
	std::random_device rd{};
	m_gen_ptr = std::unique_ptr<std::mt19937>(new std::mt19937{rd()});
	m_voltage.store(0);
	m_startSignal.store(!waitForController);
	m_state = {0,numbers::pi,0,0};
}


void RIPRealTime::runSimTime(double sim_time_step, double run_time_interval, unsigned int n_steps) {
	// start the sensor server thread
	m_log_items.reserve(n_steps/m_log_resolution + 1);
	boost::asio::io_context io_cont;
	thread sensor_server_thread(run_server_thread, std::ref(*this), std::ref(io_cont));
	array<double, 4> result_state = {0,0,0,0};
	double end_time = 0;
	auto start = chrono::steady_clock::now();
	auto nanosecond_duration = chrono::nanoseconds((int(run_time_interval * 1e9)));
	int nWait = 0;
	int nMaxWait = 800;
	while (!m_startSignal.load() && ++nWait < nMaxWait){
		this_thread::sleep_until(start + nanosecond_duration);
		start = chrono::steady_clock::now();		
	}
	if (nWait == nMaxWait){
		cerr << "control not started in time" << endl;
	}
	for (unsigned int i=0; i<n_steps; ++i){
		double voltage = m_voltage.load();
		m_rip_model.simulateWithVoltage(voltage, sim_time_step, result_state, end_time);
		{
			std::lock_guard<std::mutex> guard(m_lock);
			m_time = end_time;
			m_state = result_state;
			if (i%m_log_resolution == 0){
				ns::RotPendLogItem li;
				li.voltage = voltage;
				li.arm_angle = result_state[0];
				li.pend_angle = result_state[1];
				li.arm_ang_vel = result_state[2];
				li.pend_ang_vel = result_state[3];
				li.end_time = end_time;
				m_log_items.push_back(li);
			}
		}
		this_thread::sleep_until(start + nanosecond_duration);
		start = chrono::steady_clock::now();
	}
	cerr << "sleep at end" << endl;
	this_thread::sleep_until(start + 100*nanosecond_duration);
	cerr << "stop io" << endl;
	io_cont.stop();
	cerr << "join" << endl;
	sensor_server_thread.join();
	cerr << "exit ok" << endl;	
}

ns::RotPendSensorItem RIPRealTime::getCurrentSensorItem(){
	ns::RotPendSensorItem si;
	std::lock_guard<std::mutex> guard(m_lock);
	si.arm_angle = m_state[0];
	si.pend_angle = m_state[1];
	si.arm_angle_encoder = encoderFromAngle(si.arm_angle);
	si.pend_angle_encoder = encoderFromAngle(si.pend_angle);
	si.arm_ang_vel = m_state[2];
	si.pend_ang_vel = m_state[3];
	si.time = m_time;
	return si;
}

std::string RIPRealTime::getCurrentSensorItemJsonString(){
	json j;
	std::lock_guard<std::mutex> guard(m_lock);
	j["arm angle"] = m_state[0];
	j["pend angle"] = m_state[1];
	j["arm ang vel"] = m_state[2];
	j["pend ang vel"] = m_state[3];
	j["arm angle enc"] = encoderFromAngle(m_state[0]);
	j["pend angle enc"] = encoderFromAngle(m_state[1]);
	j["time"] = m_time;
	return j.dump();
}

void RIPRealTime::applyVoltage(double voltage){
	m_voltage.store(voltage);
}

void RIPRealTime::startSimulation(){
	m_startSignal.store(true);
}

RIPRealTime::~RIPRealTime(){
	json j;
	json log_items = m_log_items;
	j["horizontal arm length"] = m_rip_model.getHorizontalArmLength();
	j["pendulum arm length"] = m_rip_model.getPendulumArmLength();
	j["log items"] = log_items;
	std::ofstream f(m_log_file_name);
	f << setw(4) << j << std::endl;
}
