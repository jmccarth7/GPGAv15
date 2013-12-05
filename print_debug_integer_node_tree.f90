subroutine print_debug_integer_node_tree( label, input_array  )




! print INTEGER arrays of the form:

!  input_array(1:n_Nodes,1:n_Trees, 1:n_GP_individuals )


use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module

implicit none


character(*) :: label


integer(kind=4) :: i_GP_Individual

integer(kind=4) :: i_tree
integer(kind=4) :: i_node


!integer(kind=4) :: i_parm

integer(kind=4), dimension( 1:n_Nodes,1:n_Trees, 1:n_GP_individuals ) :: &
                         input_array

!--------------------------------------------------------------------------------

write(6,'(/A)') 'pd3: entry print_debug3'


!! debug
write(6,'(/A,1x,A)') 'pd3: print ', label
write(6,'(A)') &
   'pd3: i_GP_Individual, i_tree, i_node input_array(i_Node,i_Tree, i_GP_individual )'
do  i_GP_individual = 1, n_GP_individuals
    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes

            if( input_array(i_Node,i_Tree, i_GP_individual ) > -9999 )then

                write(6,'(3x,4(1x,I10))') &
                  i_GP_Individual, i_tree, i_node, &
                         input_array(i_Node,i_Tree, i_GP_individual )
            endif ! input_array(i_Node,i_Tree, i_GP_individual ) > 0.0d0

        enddo
    enddo
enddo ! i_GP_individual


write(6,'(/A)') ' '

!-------------------------------------------------------------------------------------------------

write(6,'(A//)') 'pd3: at return   '

return

end subroutine print_debug_integer_node_tree
