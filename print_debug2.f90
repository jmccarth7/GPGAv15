subroutine print_debug2( label, input_array  ) 


use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use Runge_Kutta_Variables_module

implicit none

!real(kind=8) :: cff

character(*) :: label

!integer(kind=4),dimension(n_GP_Individuals)  :: Ranked_Fitness_Index

integer(kind=4) :: i_GP_Individual
!integer(kind=4) :: j_GP_Individual

!integer(kind=4) :: icff
!integer(kind=4) :: i_tree
!integer(kind=4) :: i_node

integer(kind=4) :: i_parm

!real(kind=8), dimension( 1:n_Nodes,1:n_Trees, 1:n_GP_individuals ) :: &
!                         input_array

!real(kind=8), dimension( 1:n_CODE_equations, 1:n_GP_individuals ) :: &
!                         input_array

real(kind=8), dimension(1:n_Maximum_Number_Parameters, 1:n_GP_Individuals) :: &
                         input_array
!--------------------------------------------------------------------------------

write(6,'(/A)') 'pd2: entry print_debug1'


!! debug
!write(6,'(/A,1x,A)') 'pd2: print ', label
!write(6,'(A)') &
!   'pd2: i_GP_Individual, i_tree, i_node input_array(i_Node,i_Tree, i_GP_individual )'
!do  i_GP_individual = 1, n_GP_individuals
!    do  i_tree = 1, n_trees
!        do  i_node = 1, n_nodes
!            !if( GP_Child_Population_Node_Type(i_Node,i_Tree, i_GP_individual ) == 0 )then
!            if( input_array(i_Node,i_Tree, i_GP_individual ) > 0.0d0 )then
!                    write(6,'(1x,I6,2(1x,I3), 1x, E15.7)') &
!                      'pd2:', i_GP_Individual, i_tree, i_node, &
!                             input_array(i_Node,i_Tree, i_GP_individual )
!            endif ! GP_population_node_parameters(i_Node,i_Tree, i_GP_individual ) > 0.0d0
!            !endif ! GP_Child_Population_Node_Type(i_Node,i_Tree, i_GP_individual ) == 0
!        enddo
!    enddo
!enddo ! i_GP_individual


!write(6,'(/A)') ' '

!!! debug
write(6,'(/A,1x,A)') 'pd2: print ', label
write(6,'(A)') &
   'pd2: i_parm, input_array(i_parm, i_GP_individual )'
do  i_GP_individual = 1, n_GP_individuals
    do  i_parm = 1, n_Maximum_Number_Parameters
        if( abs( input_array(i_parm, i_GP_individual ) ) > 0.0d0 )then
            write(6,'(I6,1x,I6, 10x, E15.7)') &
                  i_GP_Individual, i_parm, &
                    input_array(i_parm, i_GP_individual )
        endif ! abs( input_array(i_parm, i_GP_individual ) ) > 0.0d0
    enddo
enddo ! i_GP_individual

!-------------------------------------------------------------------------------------------------

write(6,'(A//)') 'pd2: at return   '

return

end subroutine print_debug2
