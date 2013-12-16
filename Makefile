PROG =	GPCODE_Tree

SRCS =	0GPCODE_GA_lmdif_Parameter_Optimization_test.f90 allocate_arrays1.f90 \
	bcast1.f90 bcast2.f90 bcast3.f90 build_trees.f90 calc_fitness.f90 \
	calc_stats.f90 check_for_elite.f90 class_serialization_visitor.f90 \
	clock_module.f90 combine_tree_strings.f90 comp_data_variance.f90 \
	count_parens.f90 create_equations.f90 create_tree_node_string.f90 \
	deallocate_arrays1.f90 deserialize_trees.f90 deserialize_trees2.f90 \
	enorm.f90 Fasham_Forcing.f90 fasham_variables_module.f90 fcn.f90 \
	fdjac2.f90 fill_string_arrays.f90 \
	GA_Fitness_Proportionate_Asexual_Reproduction.f90 GA_Mutations.f90 \
	GA_parameters_module.f90 GA_random_replace.f90 \
	GA_replace_bad_individuals.f90 GA_save_elites.f90 \
	GA_Tournament_Style_Sexual_Reproduction.f90 GA_variables_module.f90 \
	Generate_Dot_Graph.f90 GP_calc_diversity_index.f90 \
	GP_calc_fitness.f90 GP_Check_Terminals.f90 GP_Clean_Tree_Nodes.f90 \
	GP_data_module.f90 GP_Fitness_Proportionate_Asexual_Reproduction.f90 \
	GP_Mutations.f90 \
	GP_para_lmdif_process.f90 GP_parameters_module.f90 \
	GP_ranking_sort.f90 GP_select_best_RK_lmdif_result.f90 \
	GP_Tournament_Style_Sexual_Reproduction.f90 GP_Tree_Build.f90 \
	GP_Tree_Build_single.f90 GP_Tree_Swap.f90 GP_variables_module.f90 \
	GPCODE_GA_lmdif_Parameter_Optimization.f90 indiv_fitness.f90 \
	init_values.f90 init_values_LV.f90 init_values_NPZ.f90 \
	Initialize_GA_Child_Parameters.f90 initialize_model.f90 lmdif.f90 \
	lmpar.f90 load_pow2_level.f90 mpi_module.f90 Numerical_methods.f90 \
	parse_fbio_strings.f90 print4.f90 print_debug_integer_node_tree.f90 \
	print_debug_real_node_tree.f90 print_debug_real_nparm.f90 \
	print_entire_tree.f90 print_trees.f90 print_values1.f90 \
	print_values2.f90 qrfac.f90 qrsolv.f90 random_real.f90 \
	read_cntl_stuff.f90 reduce_constant.f90 reduce_expression.f90 \
	remove_abs_zero.f90 remove_double_parens.f90 remove_string_blanks.f90 \
	RKBM.f90 rm_exp_paren.f90 Runge_Kutta_Box_Model_new.f90 \
	select_best_RK_lmdif_result.f90 serialize_trees.f90 \
	set_answer_arrays.f90 set_modified_indiv.f90 setup_run_fcn.f90 \
	setup_run_lmdif.f90 setup_run_para_lmdif.f90 sort.f90 sse0_calc.f90 \
	summary_GP_indiv.f90 summary_GP_indiv2.f90 swap_module.f90 \
	Tree_Helper_module.f90 class_tree_node.f03 Global_Setup.f03 \
	Interfaces.f03 Math_Node_Functions.f03  \
	tree_node_factory_module.f03

OBJS =	0GPCODE_GA_lmdif_Parameter_Optimization_test.o allocate_arrays1.o \
	bcast1.o bcast2.o bcast3.o build_trees.o calc_fitness.o calc_stats.o \
	check_for_elite.o class_serialization_visitor.o clock_module.o \
	combine_tree_strings.o comp_data_variance.o count_parens.o \
	create_equations.o create_tree_node_string.o deallocate_arrays1.o \
	deserialize_trees.o deserialize_trees2.o enorm.o Fasham_Forcing.o \
	fasham_variables_module.o fcn.o fdjac2.o fill_string_arrays.o \
	GA_Fitness_Proportionate_Asexual_Reproduction.o GA_Mutations.o \
	GA_parameters_module.o GA_random_replace.o \
	GA_replace_bad_individuals.o GA_save_elites.o \
	GA_Tournament_Style_Sexual_Reproduction.o GA_variables_module.o \
	Generate_Dot_Graph.o GP_calc_diversity_index.o GP_calc_fitness.o \
	GP_Check_Terminals.o GP_Clean_Tree_Nodes.o GP_data_module.o \
	GP_Fitness_Proportionate_Asexual_Reproduction.o \
	GP_Mutations.o GP_para_lmdif_process.o \
	GP_parameters_module.o GP_ranking_sort.o \
	GP_select_best_RK_lmdif_result.o \
	GP_Tournament_Style_Sexual_Reproduction.o GP_Tree_Build.o \
	GP_Tree_Build_single.o GP_Tree_Swap.o GP_variables_module.o \
	GPCODE_GA_lmdif_Parameter_Optimization.o indiv_fitness.o \
	init_values.o init_values_LV.o init_values_NPZ.o \
	Initialize_GA_Child_Parameters.o initialize_model.o lmdif.o lmpar.o \
	load_pow2_level.o mpi_module.o Numerical_methods.o \
	parse_fbio_strings.o print4.o print_debug_integer_node_tree.o \
	print_debug_real_node_tree.o print_debug_real_nparm.o \
	print_entire_tree.o print_trees.o print_values1.o print_values2.o \
	qrfac.o qrsolv.o random_real.o read_cntl_stuff.o reduce_constant.o \
	reduce_expression.o remove_abs_zero.o remove_double_parens.o \
	remove_string_blanks.o RKBM.o rm_exp_paren.o \
	Runge_Kutta_Box_Model_new.o select_best_RK_lmdif_result.o \
	serialize_trees.o set_answer_arrays.o set_modified_indiv.o \
	setup_run_fcn.o setup_run_lmdif.o setup_run_para_lmdif.o sort.o \
	sse0_calc.o summary_GP_indiv.o summary_GP_indiv2.o swap_module.o \
	Tree_Helper_module.o class_tree_node.o Global_Setup.o Interfaces.o \
	Math_Node_Functions.o  tree_node_factory_module.o

LIBS =	

CC = cc
CFLAGS = -O

#FC = gfortran
#FFLAGS = -g
#F90 = gfortran
#F90FLAGS = -g
#LDFLAGS = -Wl,-no_pie
#LIBS= -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk/usr/lib -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk/usr/lib  -L/Developer/SDKs/MacOSX10.6.sdk/usr/lib


# note: mpif90 is based on gfortran
FC = /opt/openmpi-1.6.5/bin/mpif90
#FFLAGS =  -O3  -ffree-form  #-fbacktrace  -fcheck=bounds  # -Wall  # -fdefault-integer-8  # -FR = -free
FFLAGS =  -g -fbacktrace -ffree-form #-fbacktrace  -fcheck=bounds  # -Wall  #-fdefault-integer-8  # -FR = -free

# note: mpif90 is based on gfortran
F90 = /opt/openmpi-1.6.5/bin/mpif90
#F90FLAGS =  -O3 -ffree-form #-fbacktrace  -fcheck=bounds  # -Wall  #-fdefault-integer-8  # -FR = -free
F90FLAGS =  -g -fbacktrace -ffree-form #-fbacktrace  -fcheck=bounds  # -Wall  #-fdefault-integer-8  # -FR = -free

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

.SUFFIXES: $(SUFFIXES) .f03

.f03.o:
	$(F90) $(F90FLAGS) -c $<

0GPCODE_GA_lmdif_Parameter_Optimization_test.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o  \
	GP_parameters_module.o GP_variables_module.o class_tree_node.o \
	mpi_module.o tree_node_factory_module.o
allocate_arrays1.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o  GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
bcast1.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o mpi_module.o
bcast2.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o mpi_module.o
bcast3.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o mpi_module.o
build_trees.o:  class_tree_node.o  tree_node_factory_module.o GP_variables_module.o Interfaces.o \
	fasham_variables_module.o  mpi_module.o 
calc_fitness.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o mpi_module.o
check_for_elite.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
class_serialization_visitor.o: GP_variables_module.o class_tree_node.o
combine_tree_strings.o:  GP_parameters_module.o \
	GP_variables_module.o
comp_data_variance.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
count_parens.o: GP_parameters_module.o GP_variables_module.o
create_equations.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
create_tree_node_string.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
deallocate_arrays1.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o  GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
deserialize_trees.o: GP_variables_module.o Tree_Helper_module.o \
	class_tree_node.o tree_node_factory_module.o
deserialize_trees2.o: GP_variables_module.o Tree_Helper_module.o \
	class_serialization_visitor.o class_tree_node.o mpi_module.o \
	tree_node_factory_module.o
Fasham_Forcing.o: GP_variables_module.o fasham_variables_module.o
fcn.o: GA_parameters_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
fill_string_arrays.o: GP_parameters_module.o GP_variables_module.o
GA_Fitness_Proportionate_Asexual_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o
GA_Mutations.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o
GA_random_replace.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
GA_replace_bad_individuals.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
GA_save_elites.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
GA_Tournament_Style_Sexual_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o
GA_variables_module.o: GA_parameters_module.o
Generate_Dot_Graph.o:  class_tree_node.o
GP_calc_diversity_index.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
GP_calc_fitness.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
GP_Check_Terminals.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
GP_Clean_Tree_Nodes.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
GP_data_module.o: GP_parameters_module.o
GP_Fitness_Proportionate_Asexual_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
GP_Mutations.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o mpi_module.o
GP_para_lmdif_process.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	clock_module.o mpi_module.o
GP_ranking_sort.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
GP_select_best_RK_lmdif_result.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o clock_module.o mpi_module.o
GP_Tournament_Style_Sexual_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
GP_Tree_Build.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
GP_Tree_Build_single.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
GP_Tree_Swap.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o mpi_module.o
GP_variables_module.o: GP_parameters_module.o \
	class_tree_node.o
GPCODE_GA_lmdif_Parameter_Optimization.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o clock_module.o mpi_module.o
indiv_fitness.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o
init_values.o: GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
init_values_LV.o: GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
init_values_NPZ.o: GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
Initialize_GA_Child_Parameters.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
initialize_model.o: GP_parameters_module.o GA_parameters_module.o GP_variables_module.o \
	fasham_variables_module.o mpi_module.o
load_pow2_level.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
mpi_module.o: 
parse_fbio_strings.o: GP_parameters_module.o GP_variables_module.o
print4.o: GA_parameters_module.o GP_parameters_module.o
print_debug_integer_node_tree.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
print_debug_real_node_tree.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
print_debug_real_nparm.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o
print_entire_tree.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
print_trees.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o
print_values1.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
print_values2.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
random_real.o: GP_parameters_module.o
read_cntl_stuff.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
reduce_constant.o: GP_parameters_module.o GP_variables_module.o
reduce_expression.o: GP_parameters_module.o GP_variables_module.o
remove_abs_zero.o: GP_parameters_module.o
remove_double_parens.o: GP_parameters_module.o
remove_string_blanks.o: GP_parameters_module.o
RKBM.o: GP_parameters_module.o GP_variables_module.o
rm_exp_paren.o: GP_parameters_module.o GP_variables_module.o
Runge_Kutta_Box_Model_new.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o Tree_Helper_module.o \
	class_serialization_visitor.o class_tree_node.o tree_node_factory_module.o \
	mpi_module.o
select_best_RK_lmdif_result.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	clock_module.o mpi_module.o
serialize_trees.o: Tree_Helper_module.o class_serialization_visitor.o \
	class_tree_node.o
set_answer_arrays.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o  GP_parameters_module.o \
	GP_variables_module.o class_tree_node.o tree_node_factory_module.o \
	mpi_module.o
set_modified_indiv.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
setup_run_fcn.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
setup_run_lmdif.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
setup_run_para_lmdif.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
sort.o: GA_parameters_module.o GP_parameters_module.o swap_module.o
sse0_calc.o: GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
summary_GP_indiv.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
summary_GP_indiv2.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	mpi_module.o
Tree_Helper_module.o: class_tree_node.o
class_tree_node.o: Math_Node_Functions.o
Global_Setup.o: Math_Node_Functions.o tree_node_factory_module.o
Interfaces.o: class_tree_node.o
tree_node_factory_module.o: class_tree_node.o
