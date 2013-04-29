subroutine deallocate_arrays1( )

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



! deallocate variable dimension arrays

deallocate( individual_elites )
deallocate( individual_SSE )
deallocate( integrated_SSE )
deallocate( individual_ranked_fitness )
deallocate( integrated_ranked_fitness )
deallocate( fitness_expectation_value )
deallocate( Run_GA_lmdif )
deallocate( Data_Array  )
deallocate( Data_Variance )

deallocate( GP_Population_Node_Parameters )
deallocate( GP_Individual_Node_Parameters )


deallocate( GP_Node_Parameters_Answer )
deallocate( GP_Node_Type_Answer )
deallocate( GP_Node_Type_for_Plotting )
deallocate( GP_Solution )
deallocate( GP_Adult_Population_Node_Type )
deallocate( GP_Child_Population_Node_Type )
deallocate( Parent_Tree_Swap_Node_Type )
deallocate( GP_Adult_Individual_SSE )
deallocate( GP_Child_Individual_SSE )
!deallocate( GP_Individual_Ranked_Fitness )
deallocate( GP_Population_Ranked_Fitness )
!deallocate( GP_Integrated_Ranked_Fitness )
deallocate( Run_GP_Calculate_Fitness )

deallocate( GA_Adult_Individual_SSE )
deallocate( GA_Child_Individual_SSE )
deallocate( GA_Integrated_SSE )
deallocate( GA_Individual_Ranked_Fitness )
deallocate( GA_Integrated_Ranked_Fitness )

deallocate( GP_Population_Parameter_Solution )
deallocate( GP_Individual_Parameter_Solution )

                                                                                                                 
deallocate( GP_Individual_N_GP_param )  ! jjm 20130409                                           
                                                                                                                 

deallocate( GP_Population_Initial_Conditions )
deallocate( GP_Individual_Initial_Conditions )

deallocate( GP_Population_Fitness )
deallocate( GP_Integrated_Population_Ranked_Fitness )

deallocate( GP_diversity_index ) 

deallocate( GP_Population_Node_Type )
deallocate( GP_Individual_Node_Type )

deallocate( Runge_Kutta_Node_Type )
deallocate( Runge_Kutta_Node_Parameters )
deallocate( Runge_Kutta_Initial_Conditions )
deallocate( Runge_Kutta_Solution )


deallocate( Node_Values )
deallocate( Tree_Evaluation )

deallocate( Tree_Value )

deallocate( Node_Eval_Type )

deallocate( bioflo )
deallocate( b_tmp )

deallocate( kval )
deallocate( btmp )
deallocate( fbio )


deallocate( node_type_string )
deallocate( node_parameters_string )
deallocate( tree_evaluation_string )
deallocate( linked_parms )


deallocate( Node_Probability )

!>>>>>>>>>>>>>                                                                                                                                   
deallocate( GP_Adult_Population_Parameter_Solution ) 
deallocate( GP_Adult_Population_SSE ) 
deallocate( GP_Child_Population_Parameter_Solution )
deallocate( GP_Child_Population_SSE )
!>>>>>>>>>>>>>                                                                                                                                   
             

return

end subroutine deallocate_arrays1
