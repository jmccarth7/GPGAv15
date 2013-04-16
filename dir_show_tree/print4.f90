subroutine print4( it, icff, &
                   left_node_value,  left_node_value_string, &
                   right_node_value, right_node_value_string, &
                   tree_evaluation )

!use parms_module
use GP_parameters_module

implicit none

real(kind=8) ::  tree_evaluation(n_nodes,n_trees)
!character(str_len), dimension( n_nodes,n_trees) ::  tree_evaluation_string

real(kind=8) ::  left_node_value,right_node_value

character(str_len) ::  left_node_value_string
character(str_len) ::  right_node_value_string

integer icff


!integer i_tree
!integer i_function

integer it

!------------------------------------------------------------------------

if( it <= n_time_steps )then
    write(6,'(8x, A, 1x,I2,3(1x,E15.7) )') &
          'icff, left, right, tree_eval ', &
           icff, left_node_value, right_node_value, &
           tree_evaluation(i_function,i_tree)
    write(6,'(8x, A, A )') 'left_node_value_string  ', &
                      trim( left_node_value_string )
    write(6,'(8x, A, A )') 'right_node_value_string ', &
                      trim( right_node_value_string )
    write(6,'(8x, A, A )') 'tree_evaluation_string  ', &
                      trim( tree_evaluation_string(i_function,i_tree) )

endif !  it <= n_time_steps

return

END subroutine print4
