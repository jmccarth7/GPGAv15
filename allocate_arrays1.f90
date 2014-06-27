subroutine allocate_arrays1( )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use kinds_mod 
use mpi
use mpi_module

use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module



implicit none

!----------------------------------------------------------------------------------------


! allocate variable dimension arrays

if( myid == 0 )then
    write(6,'(A,1x,I6)')'allo: n_code_equations = ', n_code_equations
    write(6,'(A,1x,I6)')'allo: n_nodes          = ', n_nodes
    write(6,'(A,1x,I6)')'allo: n_trees          = ', n_trees
    write(6,'(A,1x,I6)')'allo: n_levels         = ', n_levels
    write(6,'(A,1x,I6)')'allo: n_Tracked_Resources', n_Tracked_Resources
    !flush(6)
endif ! myid == 0


allocate( ga_individual_elites( n_GA_individuals )  )

allocated_memory = allocated_memory + real( n_GA_individuals * 4, kind=8 )

allocate( Run_GA_lmdif( n_GA_individuals )  ) ! orig
!!allocate( Run_GA_lmdif( 10 * n_GA_individuals )  ) ! debug only

allocated_memory = allocated_memory + real( n_GA_individuals * 4, kind=8 )

allocate( Data_Array( 0:n_time_steps, n_CODE_equations )  )

allocated_memory = allocated_memory + &
                   real( (1+n_time_steps)* n_code_equations * 8, kind=8 )


allocate( Data_Variance_inv( n_CODE_equations )  )
allocated_memory = allocated_memory + &
                   real( n_code_equations * 8, kind=8 )

allocate( ratio_Data_Variance_inv( n_CODE_equations )  )
allocated_memory = allocated_memory + &
                   real( n_code_equations * 8, kind=8 )

allocate( Parent_Tree_Swap_Node_Type(n_Nodes,2) )
allocated_memory = allocated_memory + &
                   real( n_nodes * 2 * 4     , kind=8 )

allocate( Run_GP_Calculate_Fitness(n_GP_Individuals) )
allocated_memory = allocated_memory + &
                   real( n_GP_Individuals * 4, kind=8 )


allocate( GP_Adult_Individual_SSE(n_GP_Individuals) )
allocated_memory = allocated_memory + &
                   real( n_GP_Individuals * 8, kind=8 )

allocate( GP_Child_Individual_SSE(n_GP_Individuals) )
allocated_memory = allocated_memory + &
                   real( n_GP_Individuals * 8, kind=8 )



allocate( individual_SSE( n_GA_individuals )  )
allocated_memory = allocated_memory + &
                   real( n_GA_Individuals * 8, kind=8 )


allocate( GA_Integrated_SSE(n_GA_Individuals) )
allocated_memory = allocated_memory + &
                   real( n_GA_Individuals * 8, kind=8 )

allocate( integrated_SSE( n_GA_individuals )  )
allocated_memory = allocated_memory + &
                   real( n_GA_Individuals * 8, kind=8 )


allocate( GP_n_parms( n_GP_individuals )  )
allocated_memory = allocated_memory + &
                   real( n_GP_Individuals * 4, kind=8 )


allocate( GA_Individual_Ranked_Fitness(n_GA_Individuals) )
allocated_memory = allocated_memory + &
                   real( n_GA_Individuals * 8, kind=8 )

allocate( individual_ranked_fitness( n_GA_individuals )  )
allocated_memory = allocated_memory + &
                   real( n_GA_Individuals * 8, kind=8 )


allocate( GA_Integrated_Ranked_Fitness(n_GA_Individuals) )
allocated_memory = allocated_memory + &
                   real( n_GA_Individuals * 8, kind=8 )

allocate( integrated_ranked_fitness( n_GA_individuals )  )
allocated_memory = allocated_memory + &
                   real( n_GA_Individuals * 8, kind=8 )


allocate( GP_Population_Fitness(0:n_GP_individuals) )  !00000000000000000000000000
allocated_memory = allocated_memory + &
                   real( n_GA_Individuals * 8, kind=8 )

allocate( GP_Integrated_Population_Ranked_Fitness( n_GP_Individuals ) )
allocated_memory = allocated_memory + &
                   real( n_GA_Individuals * 8, kind=8 )


allocate( GP_Individual_Ranked_Fitness(n_GP_Individuals) )
allocated_memory = allocated_memory + &
                   real( n_GP_Individuals * 8, kind=8 )

allocate( GP_Integrated_Ranked_Fitness(n_GP_Individuals) )
allocated_memory = allocated_memory + &
                   real( n_GP_Individuals * 8, kind=8 )


allocate( GP_Population_Ranked_Fitness(n_GP_Individuals) )
allocated_memory = allocated_memory + &
                   real( n_GP_Individuals * 8, kind=8 )


allocate( GP_Population_Initial_Conditions(n_CODE_equations,n_GP_individuals) )
allocated_memory = allocated_memory + &
                   real( n_code_equations * n_GP_Individuals * 8, kind=8 )


allocate( GP_Adult_Population_Node_Type( n_Nodes,n_Trees, n_GP_Individuals ) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * n_GP_Individuals * 8, kind=8 )

allocate( GP_Child_Population_Node_Type( n_Nodes,n_Trees, n_GP_Individuals ) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * n_GP_Individuals * 8, kind=8 )


allocate( GP_Population_Node_Parameters( n_nodes,n_trees, n_GP_Individuals ) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * n_GP_Individuals * 8, kind=8 )


allocate( GP_Individual_Initial_Conditions(n_CODE_equations) )
allocated_memory = allocated_memory + &
                   real( n_CODE_equations * 8, kind=8 )

allocate( GP_Individual_Node_Type(n_nodes,n_trees) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees  * 4, kind=8 )

allocate( GP_Individual_Node_Parameters(n_nodes,n_trees) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 8, kind=8 )


allocate( GP_minSSE_Individual_Initial_Conditions(n_CODE_equations) )
allocated_memory = allocated_memory + &
                   real( n_CODE_equations * 8, kind=8 )

allocate( GP_minSSE_Individual_Node_Type(n_nodes,n_trees) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 8, kind=8 )

allocate( GP_minSSE_Individual_Node_Parameters(n_nodes,n_trees) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 8, kind=8 )


allocate( GP_Individual_N_GP_param(n_GP_Individuals) )  ! jjm 20130409
allocated_memory = allocated_memory + &
                   real( n_GP_Individuals * 4, kind=8 )



allocate( GP_Node_Parameters_Answer(n_Nodes,n_Trees) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 8, kind=8 )

allocate( GP_Node_Type_Answer(n_Nodes,n_Trees) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 4, kind=8 )



allocate( GP_Node_Type_for_Plotting( n_Nodes,n_Trees, n_GP_Individuals ) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * n_GP_Individuals * 4, kind=8 )



allocate( GP_diversity_index( n_GP_individuals ) )
allocated_memory = allocated_memory + &
                   real( n_GP_Individuals * 4, kind=8 )


!---------------------------------------------------------------

allocate( Truth_Initial_Conditions( 1:n_code_equations )  )
allocated_memory = allocated_memory + &
                   real( n_CODE_equations * 8, kind=8 )

allocate( Truth_Node_Type( n_nodes, n_trees )  )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 4, kind=8 )

allocate( Truth_Node_Parameters( n_nodes, n_trees )  )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 8, kind=8 )



allocate( Truth_Model_Match( n_gp_generations ) )
allocated_memory = allocated_memory + &
                   real( n_gp_generations * 4, kind=8 )




!---------------------------------------------------------------


allocate( Node_Values(n_nodes,n_trees) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 4, kind=8 )

allocate( Tree_Evaluation(n_nodes,n_trees) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 8, kind=8 )



!allocate( GP_Trees( n_trees,  1 )  )


allocate( Tree_Value(n_trees) )
allocated_memory = allocated_memory + &
                   real( n_trees * 8, kind=8 )


allocate( Node_Eval_Type(n_nodes,n_trees) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 4, kind=8 )


allocate( Numerical_CODE_Initial_Conditions( 1:n_CODE_equations ) )
allocated_memory = allocated_memory + &
                   real( n_CODE_equations * 8, kind=8 )


allocate( Numerical_CODE_Forcing_Functions( n_CODE_forcing ) )
allocated_memory = allocated_memory + &
                   real( n_CODE_forcing * 8, kind=8 )


!allocate( Numerical_CODE_Solution( 0:n_input_data_points, n_CODE_equations ) )
allocate( Numerical_CODE_Solution( 0:n_time_steps, n_CODE_equations ) )
allocated_memory = allocated_memory + &
                   real( (1+n_time_steps) * n_CODE_equations * 8, kind=8 )



allocate( RK_Solution( 0:n_time_steps, n_CODE_equations )  )
allocated_memory = allocated_memory + &
                   real( (1+n_time_steps) * n_CODE_equations * 8, kind=8 )

allocate( RK_Node_Parameters(n_nodes,n_trees) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 8, kind=8 )

allocate( RK_Node_Type(n_nodes,n_trees) )
allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 4, kind=8 )

allocate( RK_Initial_Conditions(n_CODE_equations) )
allocated_memory = allocated_memory + &
                   real( n_CODE_equations * 8, kind=8 )


allocate( bioflo(0:n_CODE_equations,0:n_CODE_equations) )
allocated_memory = allocated_memory + &
                   real( (1+n_CODE_equations)**2 * 8, kind=8 )

allocate( bioflo_map( 1:n_CODE_equations,1:n_Tracked_Resources ) )
allocated_memory = allocated_memory + &
                   real( (1+n_CODE_equations)**2 * 4, kind=8 )


allocate( b_tmp( n_CODE_equations) )
!allocate( b_tmp(n_variables)      )
allocated_memory = allocated_memory + &
                   real( n_CODE_equations * 8, kind=8 )


allocate( GP_Trees( n_Trees, n_Tracked_Resources) )
allocated_memory = allocated_memory + &
                   real( n_Trees * n_Tracked_Resources * 8, kind=8 )


! Runge-Kutta specific work arrays

allocate( kval(4,n_CODE_equations) )
allocated_memory = allocated_memory + &
                   real( 4 * n_CODE_equations * 8, kind=8 )

allocate( btmp( n_CODE_equations) )
allocated_memory = allocated_memory + &
                   real( n_CODE_equations * 8, kind=8 )

allocate( fbio( n_CODE_equations) )
allocated_memory = allocated_memory + &
                   real( n_CODE_equations * 8, kind=8 )


!allocate( kval(4, n_variables) )
!allocate( btmp(n_variables) )
!allocate( fbio(n_variables) )

if( n_input_vars > 0 )then
    allocate( RK_data_array( 1:n_input_vars ) )
    allocated_memory = allocated_memory + &
                       real( n_input_vars * 8, kind=8 )

endif


if( L_print_equations )then

    allocate( bioflo_string(0:n_CODE_equations,0:n_CODE_equations) )
    allocated_memory = allocated_memory + &
                   real( (1+n_CODE_equations)**2 * 4, kind=8 )

    allocate( node_type_string( n_nodes, n_trees ) )
    allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 4 * 4, kind=8 )

    allocate( node_parameters_string( n_nodes, n_trees ) )
    allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * 4 * 4, kind=8 )

    allocate( tree_evaluation_string( n_nodes, n_trees ) )
    allocated_memory = allocated_memory + &
                   real( n_nodes * n_trees * str_len * 4, kind=8 )

    allocate( tree_value_string( n_trees ) )
    allocated_memory = allocated_memory + &
                   real( n_trees * str_len * 4, kind=8 )

endif ! L_print_equations


allocate( Node_Probability( n_levels ) )
allocated_memory = allocated_memory + &
               real( n_levels * 8, kind=8 )



allocate( GP_Adult_Population_SSE( n_GP_Individuals  )  )
allocated_memory = allocated_memory + &
               real( n_GP_Individuals * 8, kind=8 )


!allocate( ppex(n_Maximum_Number_Parameters,n_GA_individuals )  )


allocate( answer( n_maximum_number_parameters ) )
allocated_memory = allocated_memory + &
               real( n_maximum_number_parameters * 8, kind=8 )


allocate( output_array( n_maximum_number_parameters ) )
allocated_memory = allocated_memory + &
               real( n_maximum_number_parameters * 8, kind=8 )



!>>>>>>>>>>>>>



ga_individual_elites  = 0

Run_GA_lmdif  = .FALSE.

Data_Array  = 0.0d0
Data_Variance_inv  = 0.0d0
ratio_Data_Variance_inv  = 0.0d0

Parent_Tree_Swap_Node_Type = 0
Run_GP_Calculate_Fitness = .FALSE.


individual_SSE  = 0.0d0

GA_Integrated_SSE = 0.0d0
integrated_SSE  = 0.0d0

GP_n_parms = 0

GA_Individual_Ranked_Fitness = 0.0d0
individual_ranked_fitness  = 0.0d0

GA_Integrated_Ranked_Fitness = 0.0d0
integrated_ranked_fitness  = 0.0d0

GP_Population_Fitness = 0.0d0
GP_Integrated_Population_Ranked_Fitness = 0.0D0

GP_Individual_Ranked_Fitness = 0.0d0
GP_Integrated_Ranked_Fitness = 0.0d0
GP_Population_Ranked_Fitness = 0.0d0

GP_Population_Initial_Conditions = 0.0d0
GP_Population_Node_Parameters = 0.0d0

GP_Individual_Initial_Conditions = 0.0d0
GP_Individual_Node_Type = -9999
GP_Individual_Node_Parameters = 0.0d0

GP_minSSE_Individual_Initial_Conditions = 0.0d0
GP_minSSE_Individual_Node_Type = -9999
GP_minSSE_Individual_Node_Parameters = 0.0d0

GP_Individual_N_GP_param = 0

GP_Node_Parameters_Answer = 0.0d0
GP_Node_Type_Answer = -9999
GP_Node_Type_for_Plotting = -9999


!---------------------------------------------------------------

Truth_Initial_Conditions  = 0.0d0
Truth_Node_Type           = -9999
Truth_Node_Parameters     = 0.0d0
Truth_Model_Match         = .FALSE.

!---------------------------------------------------------------


GP_Adult_Population_Node_Type = -9999
GP_Child_Population_Node_Type = -9999

GP_Adult_Individual_SSE = 0.0d0
GP_Child_Individual_SSE = 0.0d0


Node_Values = 0.0d0
Tree_Evaluation = 0.0d0

Tree_Value = 0.0d0

Node_Eval_Type = 0

RK_Solution  = 0.0d0
RK_Node_Parameters = 0.0d0
RK_Node_Type = -9999
RK_Initial_Conditions = 0.0d0

bioflo = 0.0d0
bioflo_map = 0
b_tmp = 0.0d0
Numerical_CODE_Initial_Conditions = 0.0d0
Numerical_CODE_Forcing_Functions = 0.0d0
Numerical_CODE_Solution = 0.0d0


! Runge-Kutta specific work arrays
kval = 0.0d0
btmp = 0.0d0
fbio = 0.0d0

if( L_print_equations )then
    bioflo_string = ' '
    node_type_string = ' '
    node_parameters_string = ' '
    tree_evaluation_string = ' '
    tree_value_string = ' '
endif ! L_print_equations



Node_Probability = 0.0d0
!>>>>>>>>>>>>>
GP_Adult_Population_SSE = 0.0d0
!ppex = 0.0d0
!!>>>>>>>>>>>>>

answer       = 0.0d0
output_array = 0.0d0

if( myid == 0 )then
    write(6,'(/A,1x, E15.7/)') 'allo1: allocated_memory = ', allocated_memory
    !flush(6)
endif ! myid == 0

return

end subroutine allocate_arrays1
