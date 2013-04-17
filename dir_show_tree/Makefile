PROG =	show_tree

SRCS =	0GPCODE_GA_lmdif_Parameter_Optimization_test.f90 clock_module.f90 \
	allocate_arrays1.f90 deallocate_arrays1.f90 \
	fcn.f90 GA_parameters_module.f90 GA_variables_module.f90 \
	GP_Clean_Tree_Nodes.f90 GP_data_module.f90 \
	GP_model_parameters_module.f90 GP_parameters_module.f90 \
	GP_Tree_Build.f90 GP_Tree_Swap.f90 GP_variables_module.f90 \
	init_values.f90 init_values_LV.f90 init_values_NPZ.f90 \
	Initialize_GA_Child_Parameters.f90  print4.f90 \
	random_real.f90 random_real2.f90 read_cntl_stuff.f90 \
	Runge_Kutta_Box_Model.f90 Runge_Kutta_Variables_module.f90 \
	show_tree.f90 sort.f90  \
	swap_module.f90

OBJS =	0GPCODE_GA_lmdif_Parameter_Optimization_test.o clock_module.o fcn.o \
	allocate_arrays1.o deallocate_arrays1.o \
	GA_parameters_module.o GA_variables_module.o GP_Clean_Tree_Nodes.o \
	GP_data_module.o GP_model_parameters_module.o GP_parameters_module.o \
	GP_Tree_Build.o GP_Tree_Swap.o GP_variables_module.o init_values.o \
	init_values_LV.o init_values_NPZ.o Initialize_GA_Child_Parameters.o \
	 print4.o random_real.o random_real2.o read_cntl_stuff.o \
	Runge_Kutta_Box_Model.o Runge_Kutta_Variables_module.o \
	show_tree.o sort.o  swap_module.o

LIBS =	

CC = cc
CFLAGS = -O
FC = g95
FFLAGS = -g -ftrace=full -Wall
F90 = g95
F90FLAGS = -g -ftrace=full -Wall
LDFLAGS = -Wl,-no_pie
LIBS= -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk/usr/lib -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk/usr/lib  -L/Developer/SDKs/MacOSX10.6.sdk/usr/lib


all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

new: $(PROG)

	rm -f  $(OBJS) *.mod

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

0GPCODE_GA_lmdif_Parameter_Optimization_test.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o Runge_Kutta_Variables_module.o
allocate_arrays1.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o Runge_Kutta_Variables_module.o
deallocate_arrays1.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o Runge_Kutta_Variables_module.o
fcn.o: GP_data_module.o GP_parameters_module.o GP_variables_module.o
GA_variables_module.o: GA_parameters_module.o
GP_Clean_Tree_Nodes.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
GP_data_module.o: GP_parameters_module.o
GP_parameters_module.o: GP_model_parameters_module.o
GP_Tree_Build.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
GP_Tree_Swap.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
GP_variables_module.o: GP_parameters_module.o
init_values.o: GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o
init_values_LV.o: GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o
init_values_NPZ.o: GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o
Initialize_GA_Child_Parameters.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o
print4.o: GP_parameters_module.o
read_cntl_stuff.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
Runge_Kutta_Box_Model.o: GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o
Runge_Kutta_Variables_module.o: GP_parameters_module.o
show_tree.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o
sort.o: swap_module.o
