#######################################################
#
#   inverted pendulum control
#
#######################################################

CC = g++ 
FLAGS = -Wall -pedantic -Werror -Wno-psabi -std=c++20 -pthread -I/usr/include/eigen3 -O2

SOURCES = rip_controller.cpp\
	rip_kf_estimator_sr_quanser.cpp\
	rot_pend_sensor_item.cpp\
	dl_syscalls.c


MAIN_DEADLINE_BW = rip_control_main_deadline_bw.cpp

HEADERS = rip_controller.h\
	rip_kf_estimator_sr_quanser.h\
	rot_pend_sensor_item.h\
	dl_syscalls.h


CONTROLLER_FILE_DEADLINE_BW = ripControlDeadlineBW

MAKE = make

.PHONY : 

all :  $(CONTROLLER_FILE_DEADLINE_BW)

debug : FLAGS += -DDEBUG -g
debug: $(CONTROLLER_FILE_DEADLINE_BW)


clean :
	rm -f $(CONTROLLER_FILE)
	rm -f $(CONTROLLER_FILE_DEADLINE)
	rm -f $(CONTROLLER_FILE_DEADLINE_BW)
	rm -rf *.o	


$(CONTROLLER_FILE_DEADLINE_BW) : $(SOURCES) $(HEADERS) $(MAIN_DEADLINE_BW)
	$(CC) $(SOURCES) $(MAIN_DEADLINE_BW) $(FLAGS) -o $(CONTROLLER_FILE_DEADLINE_BW)

