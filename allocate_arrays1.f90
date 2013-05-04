subroutine allocate_arrays1( )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use Runge_Kutta_Variables_module


implicit none


!logical bad_node


!integer :: i
!integer :: nop
!integer :: message_len

!logical :: Lplot


!real(kind=8), allocatable, dimension(:) :: answer
!real(kind=8), allocatable, dimension(:) :: output_array

!real (kind=8) :: ssum, ssum2, totobs, dff

!real (kind=8) :: xcount

!----------------------------------------------------------------------------------------


! allocate variable dimension arrays


allocate( individual_elites( n_GA_individuals )  )
allocate( fitness_expectation_value( n_GA_individuals )  )

allocate( Run_GA_lmdif( n_GA_individuals )  )

allocate( Data_Array( 0:n_time_steps, n_CODE_equations )  )
allocate( Data_Variance( n_CODE_equations )  )

allocate( Parent_Tree_Swap_Node_Type(n_Nodes,2) )
allocate( Run_GP_Calculate_Fitness(n_GP_Individuals) )

allocate( GA_Adult_Individual_SSE(n_GA_Individuals) )
allocate( GA_Child_Individual_SSE(n_GA_Individuals) )

allocate( GP_Adult_Individual_SSE(n_GP_Individuals) )
allocate( GP_Child_Individual_SSE(n_GP_Individuals) )

allocate( individual_SSE( n_GA_individuals )  )

allocate( GA_Integrated_SSE(n_GA_Individuals) )
allocate( integrated_SSE( n_GA_individuals )  )

allocate( GA_Individual_Ranked_Fitness(n_GA_Individuals) )
allocate( individual_ranked_fitness( n_GA_individuals )  )

allocate( GA_Integrated_Ranked_Fitness(n_GA_Individuals) )
allocate( integrated_ranked_fitness( n_GA_individuals )  )


allocate( GP_Population_Parameter_Solution(n_GP_individuals,n_maximum_number_parameters) )
allocate( GP_Population_Initial_Conditions(n_GP_individuals,n_CODE_equations) )

allocate( GP_Population_Node_Type(n_GP_individuals,n_nodes,n_trees) )
allocate( GP_Adult_Population_Node_Type(n_GP_Individuals,n_Nodes,n_Trees) )
allocate( GP_Child_Population_Node_Type(n_GP_Individuals,n_Nodes,n_Trees) )

allocate( GP_Population_Node_Parameters(n_GP_individuals,n_nodes,n_trees) )
allocate( GP_Population_Fitness(n_GP_individuals) )
allocate( GP_Integrated_Population_Ranked_Fitness( n_GP_Individuals ) )

allocate( GP_Individual_Parameter_Solution(n_maximum_number_parameters) )
allocate( GP_Individual_Initial_Conditions(n_CODE_equations) )
allocate( GP_Individual_Node_Type(n_nodes,n_trees) )
allocate( GP_Individual_Node_Parameters(n_nodes,n_trees) )


allocate( GP_Individual_N_GP_param(n_GP_Individuals) )  ! jjm 20130409


allocate( GP_Node_Parameters_Answer(n_Nodes,n_Trees) )
allocate( GP_Node_Type_Answer(n_Nodes,n_Trees) )

!allocate( GP_Node_Type_for_Plotting(9, n_Nodes,n_Trees) )
allocate( GP_Node_Type_for_Plotting(n_GP_individuals, n_Nodes,n_Trees) )

allocate( GP_Solution(0:n_Time_Steps,n_CODE_Equations) )

allocate( GP_diversity_index( n_GP_individuals ) )



allocate( GP_Individual_Ranked_Fitness(n_GP_Individuals) )
allocate( GP_Population_Ranked_Fitness(n_GP_Individuals) )
allocate( GP_Integrated_Ranked_Fitness(n_GP_Individuals) )

allocate( Node_Values(n_nodes,n_trees) )
allocate( Tree_Evaluation(n_nodes,n_trees) )

allocate( Tree_Value(n_trees) )

allocate( Node_Eval_Type(n_nodes,n_trees) )

allocate( Runge_Kutta_Solution( 0:n_time_steps, n_CODE_equations )  )
allocate( Runge_Kutta_Node_Parameters(n_nodes,n_trees) )
allocate( Runge_Kutta_Node_Type(n_nodes,n_trees) )
allocate( Runge_Kutta_Initial_Conditions(n_CODE_equations) )

allocate( bioflo(0:n_CODE_equations,0:n_CODE_equations) )
allocate( b_tmp(n_CODE_equations) )

! Runge-Kutta specific work arrays
allocate( kval(4,n_CODE_equations) )
allocate( btmp(n_CODE_equations) )
allocate( fbio(n_CODE_equations) )


allocate( node_type_string( n_nodes, n_trees ) )
allocate( node_parameters_string( n_nodes, n_trees ) )
allocate( tree_evaluation_string( n_nodes, n_trees ) )

!allocate( answer( n_maximum_number_parameters ) )
!allocate( output_array( n_maximum_number_parameters ) )

allocate( linked_parms( 2, n_linked_parms_dim ) )

allocate( Node_Probability( n_levels ) )

!>>>>>>>>>>>>>
allocate( GP_Adult_Population_Parameter_Solution( 1:n_GP_Individuals,1:n_Maximum_Number_Parameters) )
allocate( GP_Child_Population_Parameter_Solution( 1:n_GP_Individuals,1:n_Maximum_Number_Parameters ) )
allocate( GP_Adult_Population_SSE( n_GP_Individuals  )  )
allocate( GP_Child_Population_SSE( n_GP_Individuals  )  )
!>>>>>>>>>>>>>



individual_elites  = 0
fitness_expectation_value  = 0.0d0

Run_GA_lmdif  = .FALSE.

Data_Array  = 0.0d0
Data_Variance  = 0.0d0

Parent_Tree_Swap_Node_Type = 0
Run_GP_Calculate_Fitness = .FALSE.

GA_Adult_Individual_SSE = 0.0d0
GA_Child_Individual_SSE = 0.0d0

individual_SSE  = 0.0d0

GA_Integrated_SSE = 0.0d0
integrated_SSE  = 0.0d0

GA_Individual_Ranked_Fitness = 0.0d0
individual_ranked_fitness  = 0.0d0

GA_Integrated_Ranked_Fitness = 0.0d0
integrated_ranked_fitness  = 0.0d0


GP_Population_Parameter_Solution = 0.0d0
GP_Population_Initial_Conditions = 0.0d0
GP_Population_Node_Type = -9999
GP_Population_Node_Parameters = 0.0d0
GP_Population_Fitness = 0.0d0
GP_Integrated_Population_Ranked_Fitness = 0.0D0

GP_Individual_Parameter_Solution = 0.0d0
GP_Individual_Initial_Conditions = 0.0d0
GP_Individual_Node_Type = -9999
GP_Individual_Node_Parameters = 0.0d0

GP_Individual_N_GP_param = 0

GP_Node_Parameters_Answer = 0.0d0
GP_Node_Type_Answer = -9999
GP_Node_Type_for_Plotting = -9999

GP_Solution = 0.0d0

GP_Adult_Population_Node_Type = -9999
GP_Child_Population_Node_Type = -9999

GP_Adult_Individual_SSE = 0.0d0
GP_Child_Individual_SSE = 0.0d0

GP_Individual_Ranked_Fitness = 0.0d0
GP_Population_Ranked_Fitness = 0.0d0
GP_Integrated_Ranked_Fitness = 0.0d0

Node_Values = 0.0d0
Tree_Evaluation = 0.0d0

Tree_Value = 0.0d0

Node_Eval_Type = 0

Runge_Kutta_Solution  = 0.0d0
Runge_Kutta_Node_Parameters = 0.0d0
Runge_Kutta_Node_Type = -9999
Runge_Kutta_Initial_Conditions = 0.0d0

bioflo = 0.0d0
b_tmp = 0.0d0

! Runge-Kutta specific work arrays
kval = 0.0d0
btmp = 0.0d0
fbio = 0.0d0


node_type_string = ' '
node_parameters_string = ' '
tree_evaluation_string = ' '


linked_parms = 0

Node_Probability = 0.0d0
!>>>>>>>>>>>>>
GP_Adult_Population_Parameter_Solution = 0.0d0
GP_Adult_Population_SSE = 0.0d0
GP_Child_Population_Parameter_Solution = 0.0d0
GP_Child_Population_SSE = 0.0d0
!>>>>>>>>>>>>>


return

end subroutine allocate_arrays1
