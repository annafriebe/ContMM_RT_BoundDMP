# Edited by Anna Friebe, building only cli_solver
UNAME := $(shell uname)

CC = g++
CXXFLAGS = -Wall -Wextra -pedantic -g -std=c++11 -O3

ifeq ($(UNAME), Linux)
CPPFLAGS = -I/usr/include/eigen3 -I./utils/tinyxml2
endif

ifeq ($(UNAME), Darwin)
CPPFLAGS = -I/opt/local/include/eigen3 -I./utils/tinyxml2
endif

LDFLAGS = -L./ 
LDLIBS = -lprosit

APP1 = cli_solver
APPS = $(APP1) 
PROSIT_OBJS = ./utils/auxiliary_functions.o \
              ./utils/help_messages.o \
              ./utils/distribution.o \
              ./utils/pmf.o \
              ./utils/cdf.o \
              ./tasks/generic_task_descriptor.o \
              ./tasks/rr_task_descriptor.o \
              ./tasks/fp_task_descriptor.o \
              ./tasks/fp_task_schedule.o \
              ./solver/probability_solver.o \
              ./solver/rr_probability_solver.o \
              ./solver/fp_probability_solver.o \
              ./solver/analytic_rr_probability_solver.o \
              ./solver/companion_rr_probability_solver.o \
              ./solver/qbd_rr_probability_solver.o \
              ./solver/lt_qbd_rr_probability_solver.o \
              ./solver/cr_qbd_rr_probability_solver.o \
              ./qos/linear_qos_function.o \
              ./qos/quadratic_qos_function.o 

CLI_SOLVER_OBJS = ./cli/cli_solver.o

all: $(APPS)

libprosit.a: $(PROSIT_OBJS)
	ar rc $@ $(PROSIT_OBJS)
	ranlib $@

$(APP1): $(CLI_SOLVER_OBJS) libprosit.a
	$(CC) ${CXXFLAGS} $(CLI_SOLVER_OBJS) -o cli_solver $(LDLIBS) $(LDFLAGS)

.PHONY: clean

clean: 
	rm -f ./solver/*.o ./solver/*.d ./cli/*.o ./cli/*.d ./qos/*.o ./qos/*.d ./tasks/*.o ./tasks/*.d ./xml/*.o ./xml/*.d ./utils/*.o ./utils/*.d ./utils/tinyxml2/*.o ./utils/tinyxml2/*.d *.a *~ $(APPS)

%.d: %.cpp
	$(CC) $(CXXFLAGS) -MM -MF $@ $<

-include $(PROSIT_OBJS:.o=.d)
-include $(CLI_SOLVER_OBJS:.o=.d)
