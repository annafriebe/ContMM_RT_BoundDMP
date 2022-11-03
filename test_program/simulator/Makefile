#######################################################
#
#   Rotary inverted pendulum simulation
#
#######################################################

CC = g++ 
FLAGS = -Wall -pedantic -Werror -Wno-psabi -std=c++20 -pthread -latomic -O2

SOURCES = quanser_rotary_inverted_pendulum_model.cpp\
	rip_real_time.cpp

MAIN_REAL_TIME= rip_real_time_main.cpp

HEADERS = quanser_rotary_inverted_pendulum_model.h\
	rip_real_time.h

REAL_TIME_FILE = ripRTSimulate

MAKE = make

.PHONY : 

all : $(REAL_TIME_FILE) $(REAL_TIME_FILE_DEADLINE)

debug : FLAGS += -DDEBUG -g
debug: $(REAL_TIME_FILE) $(REAL_TIME_FILE_DEADLINE)

 
clean :
	rm -f $(REAL_TIME_FILE)
	rm -f $(REAL_TIME_FILE_DEADLINE)
	rm -rf *.o	

$(REAL_TIME_FILE) : $(SOURCES) $(HEADERS) $(MAIN_REAL_TIME)
	$(CC) $(SOURCES) $(MAIN_REAL_TIME) $(FLAGS) -o $(REAL_TIME_FILE)

