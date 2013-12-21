subroutine allocate_arrays1( )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module

!!!!!!!!!use GP_model_parameters_module
use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module



implicit none

!----------------------------------------------------------------------------------------


! allocate variable dimension arrays


allocate( ga_individual_elites( n_GA_individuals )  )
!allocate( fitness_expectation_value( n_GA_individuals )  )

allocate( Run_GA_lmdif( n_GA_individuals )  )

allocate( Data_Array( 0:n_time_steps, n_CODE_equations )  )
allocate( Data_Variance( n_CODE_equations )  )
allocate( Data_Variance_inv( n_CODE_equations )  )

allocate( Parent_Tree_Swap_Node_Type(n_Nodes,2) )
allocate( Run_GP_Calculate_Fitness(n_GP_Individuals) )


allocate( GP_Adult_Individual_SSE(n_GP_Individuals) )
allocate( GP_Child_Individual_SSE(n_GP_Individuals) )


allocate( individual_SSE( n_GA_individuals )  )

allocate( GA_Integrated_SSE(n_GA_Individuals) )
allocate( integrated_SSE( n_GA_individuals )  )

allocate( GP_n_parms( n_GP_individuals )  )

allocate( GA_Individual_Ranked_Fitness(n_GA_Individuals) )
allocate( individual_ranked_fitness( n_GA_individuals )  )

allocate( GA_Integrated_Ranked_Fitness(n_GA_Individuals) )
allocate( integrated_ranked_fitness( n_GA_individuals )  )

allocate( GP_Population_Fitness(n_GP_individuals) )
allocate( GP_Integrated_Population_Ranked_Fitness( n_GP_Individuals ) )

allocate( GP_Individual_Ranked_Fitness(n_GP_Individuals) )
allocate( GP_Integrated_Ranked_Fitness(n_GP_Individuals) )

allocate( GP_Population_Ranked_Fitness(n_GP_Individuals) )

allocate( GP_Population_Initial_Conditions(n_CODE_equations,n_GP_individuals) )

allocate( GP_Adult_Population_Node_Type( n_Nodes,n_Trees, n_GP_Individuals ) )
allocate( GP_Child_Population_Node_Type( n_Nodes,n_Trees, n_GP_Individuals ) )

allocate( GP_Population_Node_Parameters( n_nodes,n_trees, n_GP_Individuals ) )

allocate( GP_Individual_Initial_Conditions(n_CODE_equations) )
allocate( GP_Individual_Node_Type(n_nodes,n_trees) )
allocate( GP_Individual_Node_Parameters(n_nodes,n_trees) )

allocate( child_node_parameters( n_nodes, n_trees, n_GP_Individuals ) )

allocate( GP_Individual_N_GP_param(n_GP_Individuals) )  ! jjm 20130409


allocate( GP_Node_Parameters_Answer(n_Nodes,n_Trees) )
allocate( GP_Node_Type_Answer(n_Nodes,n_Trees) )

!allocate( GP_Node_Type_for_Plotting(9, n_Nodes,n_Trees) )

allocate( GP_Node_Type_for_Plotting( n_Nodes,n_Trees, n_GP_Individuals ) )

!allocate( GP_Solution(0:n_Time_Steps,n_CODE_Equations) )

allocate( GP_diversity_index( n_GP_individuals ) )




allocate( Node_Values(n_nodes,n_trees) )
allocate( Tree_Evaluation(n_nodes,n_trees) )

allocate( Tree_Value(n_trees) )

allocate( Node_Eval_Type(n_nodes,n_trees) )

allocate( Runge_Kutta_Solution( 0:n_time_steps, n_CODE_equations )  )
allocate( Runge_Kutta_Node_Parameters(n_nodes,n_trees) )
allocate( Runge_Kutta_Node_Type(n_nodes,n_trees) )
allocate( Runge_Kutta_Initial_Conditions(n_CODE_equations) )

allocate( RK_Solution( 0:n_time_steps, n_CODE_equations )  )
allocate( RK_Node_Parameters(n_nodes,n_trees) )
allocate( RK_Node_Type(n_nodes,n_trees) )
allocate( RK_Initial_Conditions(n_CODE_equations) )

allocate( bioflo(0:n_CODE_equations,0:n_CODE_equations) )
allocate( bioflo_map( 1:n_CODE_equations,1:n_Tracked_Resources ) )
allocate( bioflo_string(0:n_CODE_equations,0:n_CODE_equations) )

!allocate( b_tmp(n_CODE_equations) )
allocate( b_tmp(n_variables)      )

allocate( Numerical_CODE_Initial_Conditions( 1:n_CODE_equations ) ) 

allocate( Numerical_CODE_Forcing_Functions( n_CODE_forcing ) )

allocate( Numerical_CODE_Solution( 0:n_time_steps, n_CODE_equations ) )
!!allocate( GP_Trees( n_Trees, n_Tracked_Resources) ) 

! Runge-Kutta specific work arrays
!allocate( kval(4,n_CODE_equations) )
!allocate( btmp(n_CODE_equations) )
!allocate( fbio(n_CODE_equations) )

allocate( kval(4, n_variables) )
allocate( btmp(n_variables) )
allocate( fbio(n_variables) )



allocate( node_type_string( n_nodes, n_trees ) )
allocate( node_parameters_string( n_nodes, n_trees ) )
allocate( tree_evaluation_string( n_nodes, n_trees ) )
allocate( tree_value_string( n_trees ) )

!!allocate( linked_parms( 2, n_linked_parms_dim ) )

allocate( Node_Probability( n_levels ) )

!>>>>>>>>>>>>>

allocate( GP_Adult_Population_SSE( n_GP_Individuals  )  )
!allocate( GP_Child_Population_SSE( n_GP_Individuals  )  )

allocate( ppe(n_Maximum_Number_Parameters,n_GA_individuals )  )


!allocate( buffer( n_maximum_number_parameters + 2 ) )
!allocate( buffer_recv( n_maximum_number_parameters + 2 ) )

!allocate( buffer2( n_maximum_number_parameters + 2 ) )
!allocate( buffer2_recv( n_maximum_number_parameters + 2 ) )



!>>>>>>>>>>>>>



ga_individual_elites  = 0
!fitness_expectation_value  = 0.0d0

Run_GA_lmdif  = .FALSE.

Data_Array  = 0.0d0
Data_Variance  = 0.0d0
Data_Variance_inv  = 0.0d0

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

GP_Individual_N_GP_param = 0

GP_Node_Parameters_Answer = 0.0d0
GP_Node_Type_Answer = -9999
GP_Node_Type_for_Plotting = -9999

!GP_Solution = 0.0d0

GP_Adult_Population_Node_Type = -9999
GP_Child_Population_Node_Type = -9999

GP_Adult_Individual_SSE = 0.0d0
GP_Child_Individual_SSE = 0.0d0


Node_Values = 0.0d0
Tree_Evaluation = 0.0d0

Tree_Value = 0.0d0

Node_Eval_Type = 0

Runge_Kutta_Solution  = 0.0d0
Runge_Kutta_Node_Parameters = 0.0d0
Runge_Kutta_Node_Type = -9999
Runge_Kutta_Initial_Conditions = 0.0d0

RK_Solution  = 0.0d0
RK_Node_Parameters = 0.0d0
RK_Node_Type = -9999
RK_Initial_Conditions = 0.0d0

bioflo = 0.0d0
bioflo_map = 0.0d0
b_tmp = 0.0d0
Numerical_CODE_Initial_Conditions = 0.0d0
Numerical_CODE_Forcing_Functions = 0.0d0
Numerical_CODE_Solution = 0.0d0

! GP_Trees = ?

! Runge-Kutta specific work arrays
kval = 0.0d0
btmp = 0.0d0
fbio = 0.0d0


node_type_string = ' '
node_parameters_string = ' '
tree_evaluation_string = ' '
tree_value_string = ' '

!!linked_parms = 0

Node_Probability = 0.0d0
!>>>>>>>>>>>>>
GP_Adult_Population_SSE = 0.0d0
!GP_Child_Population_SSE = 0.0d0
!>>>>>>>>>>>>>

!buffer       = 0.0d0
!buffer_recv  = 0.0d0
!buffer2      = 0.0d0
!buffer2_recv = 0.0d0

return

end subroutine allocate_arrays1
