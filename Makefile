PROG =	GPGACODE_test

SRCS =	0GPCODE_GA_lmdif_Parameter_Optimization_test.f90 allocate_arrays1.f90 \
	calc_fitness.f90 check_for_elite.f90 clock_module.f90 \
	comp_data_variance.f90 comp_GP_child_indiv_sse.f90 \
	deallocate_arrays1.f90 enorm.f90 fcn.f90 fdjac2.f90 \
	GA_Fitness_Proportionate_Asexual_Reproduction.f90 GA_Mutations.f90 \
	GA_parameters_module.f90 GA_replace_bad_individuals.f90 \
	GA_save_elites.f90 GA_Tournament_Style_Sexual_Reproduction.f90 \
	GA_variables_module.f90 gaussian_random_number_generator.f90 \
	GP_calc_diversity_index.f90 GP_calc_fitness.f90 \
	GP_Check_Terminals.f90 GP_Clean_Tree_Nodes.f90 GP_data_module.f90 \
	GP_Elitists.f90  \
	GP_Fitness_Proportionate_Asexual_Reproduction.f90 \
	GP_model_parameters_module.f90 GP_Mutations.f90 \
	GP_parameters_module.f90 GP_Tournament_Style_Sexual_Reproduction.f90 \
	GP_Tree_Build.f90 GP_Tree_Swap.f90 GP_variables_module.f90 \
	GPCODE_GA_lmdif_Parameter_Optimization.f90 indiv_fitness.f90 \
	init_values.f90 init_values_LV.f90 init_values_NPZ.f90 \
	Initialize_GA_Child_Parameters.f90 lmdif.f90 lmpar.f90 \
	median_calc.f90 mpi_module.f90 print4.f90 print_trees.f90 qrfac.f90 \
	qrsolv.f90 random_real.f90 read_cntl_stuff.f90 \
	Runge_Kutta_Box_Model.f90 Runge_Kutta_Variables_module.f90 \
	set_answer_arrays.f90 setup_run_fcn.f90 setup_run_lmdif.f90 sort.f90 \
	sse0_calc.f90 summary_GP_indiv.f90 swap_module.f90

OBJS =	0GPCODE_GA_lmdif_Parameter_Optimization_test.o allocate_arrays1.o \
	calc_fitness.o check_for_elite.o clock_module.o comp_data_variance.o \
	comp_GP_child_indiv_sse.o deallocate_arrays1.o enorm.o fcn.o fdjac2.o \
	GA_Fitness_Proportionate_Asexual_Reproduction.o GA_Mutations.o \
	GA_parameters_module.o GA_replace_bad_individuals.o GA_save_elites.o \
	GA_Tournament_Style_Sexual_Reproduction.o GA_variables_module.o \
	gaussian_random_number_generator.o GP_calc_diversity_index.o \
	GP_calc_fitness.o GP_Check_Terminals.o GP_Clean_Tree_Nodes.o \
	GP_data_module.o GP_Elitists.o  \
	GP_Fitness_Proportionate_Asexual_Reproduction.o \
	GP_model_parameters_module.o GP_Mutations.o GP_parameters_module.o \
	GP_Tournament_Style_Sexual_Reproduction.o GP_Tree_Build.o \
	GP_Tree_Swap.o GP_variables_module.o \
	GPCODE_GA_lmdif_Parameter_Optimization.o indiv_fitness.o \
	init_values.o init_values_LV.o init_values_NPZ.o \
	Initialize_GA_Child_Parameters.o lmdif.o lmpar.o median_calc.o \
	mpi_module.o print4.o print_trees.o qrfac.o qrsolv.o random_real.o \
	read_cntl_stuff.o Runge_Kutta_Box_Model.o \
	Runge_Kutta_Variables_module.o set_answer_arrays.o setup_run_fcn.o \
	setup_run_lmdif.o sort.o sse0_calc.o summary_GP_indiv.o swap_module.o

LIBS =	

CC = cc
CFLAGS = -O
#FC = g95
#FFLAGS = -g
#F90 = g95
#F90FLAGS = -g
#LDFLAGS = -lSystemStubs

# note: mpif90 is based on gfortran
FC = /opt/openmpi-1.6/bin/mpif90
FFLAGS =   -g -Wall  -ffree-form # -fbacktrace # -fdefault-integer-8  # -FR = -free

# note: mpif90 is based on gfortran
F90 = /opt/openmpi-1.6/bin/mpif90
F90FLAGS =   -g -Wall  -ffree-form #  -fbacktrace #-fdefault-integer-8  # -FR = -free

LDFLAGS = -L/opt/openmpi-1.6/lib \
          -I/Developer/SDKs/MacOSX10.6.sdk/usr/include
LIBS= -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk/usr/lib \
      -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk/usr/lib \
      -L/Developer/SDKs/MacOSX10.6.sdk/usr/lib



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
	GP_variables_module.o Runge_Kutta_Variables_module.o mpi_module.o
allocate_arrays1.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o mpi_module.o
calc_fitness.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o mpi_module.o
check_for_elite.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
comp_data_variance.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o mpi_module.o
comp_GP_child_indiv_sse.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o mpi_module.o
deallocate_arrays1.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o mpi_module.o
fcn.o: GA_parameters_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
GA_Fitness_Proportionate_Asexual_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o
GA_Mutations.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o
GA_replace_bad_individuals.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
GA_save_elites.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
GA_Tournament_Style_Sexual_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o
GA_variables_module.o: GA_parameters_module.o
GP_calc_diversity_index.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
GP_calc_fitness.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o mpi_module.o
GP_Check_Terminals.o: GA_parameters_module.o GA_variables_module.o \
	GP_model_parameters_module.o GP_parameters_module.o \
	GP_variables_module.o
GP_Clean_Tree_Nodes.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
GP_data_module.o: GP_parameters_module.o
GP_Elitists.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
GP_Fitness_Proportionate_Asexual_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_parameters_module.o GP_variables_module.o
GP_Mutations.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
GP_parameters_module.o: GP_model_parameters_module.o
GP_Tournament_Style_Sexual_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_parameters_module.o GP_variables_module.o
GP_Tree_Build.o: GA_parameters_module.o GA_variables_module.o \
	GP_model_parameters_module.o GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
GP_Tree_Swap.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
GP_variables_module.o: GP_parameters_module.o
GPCODE_GA_lmdif_Parameter_Optimization.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o clock_module.o mpi_module.o
indiv_fitness.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
init_values.o: GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o mpi_module.o
init_values_LV.o: GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o mpi_module.o
init_values_NPZ.o: GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o mpi_module.o
Initialize_GA_Child_Parameters.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
median_calc.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o mpi_module.o
mpi_module.o: 
print4.o: GA_parameters_module.o GP_parameters_module.o
print_trees.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o mpi_module.o
random_real.o: GP_parameters_module.o
read_cntl_stuff.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
Runge_Kutta_Box_Model.o: GA_parameters_module.o GP_parameters_module.o \
	GP_variables_module.o Runge_Kutta_Variables_module.o mpi_module.o
Runge_Kutta_Variables_module.o: GP_parameters_module.o
set_answer_arrays.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o mpi_module.o
setup_run_fcn.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
setup_run_lmdif.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
sort.o: GA_parameters_module.o GP_parameters_module.o swap_module.o
sse0_calc.o: GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
summary_GP_indiv.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	Runge_Kutta_Variables_module.o mpi_module.o
