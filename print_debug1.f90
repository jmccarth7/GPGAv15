subroutine print_debug1( label, input_array ) 



! print REAL arrays of the form:
!   input_array( 1:n_Nodes,1:n_Trees, 1:n_GP_individuals )


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
integer(kind=4) :: i_tree
integer(kind=4) :: i_node


!integer(kind=4) :: i_parm

real(kind=8), dimension( 1:n_Nodes,1:n_Trees, 1:n_GP_individuals ) :: &
                         input_array

!real(kind=8), dimension( 1:n_CODE_equations, 1:n_GP_individuals ) :: &
!                         GP_Population_Initial_Conditions_temp

!real(kind=8), dimension(1:n_Maximum_Number_Parameters, 1:n_GP_Individuals) :: &
!                         GP_Adult_Population_Parameter_Solution_temp
!--------------------------------------------------------------------------------

write(6,'(/A)') 'pd1: entry print_debug1'


!! debug
write(6,'(/A,1x,A)') 'pd1: print ', label
write(6,'(A)') &
   'pd1: i_GP_Individual, i_tree, i_node input_array(i_Node,i_Tree, i_GP_individual )'
do  i_GP_individual = 1, n_GP_individuals
    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes

            if( input_array(i_Node,i_Tree, i_GP_individual ) > 0.0d0 )then
                    write(6,'(3x,3(1x,I10), 1x, E15.7)') &
                      i_GP_Individual, i_tree, i_node, &
                             input_array(i_Node,i_Tree, i_GP_individual )
            endif ! input_array(i_Node,i_Tree, i_GP_individual ) > 0.0d0

        enddo
    enddo
enddo ! i_GP_individual


write(6,'(/A)') ' '

!!! debug
!write(6,'(/A)') 'pd1: before applying  sort to GP_Adult_Population_Parameter_Solution '
!do  i_GP_individual = 1, n_GP_individuals
!    do  i_parm = 1, n_Maximum_Number_Parameters
!        if( abs( GP_Adult_Population_Parameter_Solution(i_parm, i_GP_individual ) ) > 0.0d0 )then
!            write(6,'(A,1x,I6,1x,I6, 10x, E15.7)') &
!              'pd1: before:', i_GP_Individual, i_parm, &
!                    GP_Adult_Population_Parameter_Solution(i_parm, i_GP_individual )
!        endif ! abs( GP_Adult_Population_Parameter_Solution(i_parm, i_GP_individual ) ) > 0.0d0
!    enddo
!enddo ! i_GP_individual


!-------------------------------------------------------------------------------------------------

write(6,'(A//)') 'pd1: at return   '

return

end subroutine print_debug1
