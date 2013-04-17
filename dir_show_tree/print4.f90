subroutine print4( it, icff, &
                   left_node_value,  left_node_value_string, &
                   right_node_value, right_node_value_string, &
                   tree_evaluation )

!use parms_module
use GP_parameters_module

implicit none

real(kind=8) ::  tree_evaluation(n_nodes,n_trees)

real(kind=8) ::  left_node_value
real(kind=8) ::  right_node_value

character(str_len) ::  left_node_value_string
character(str_len) ::  right_node_value_string

integer icff


integer it

!------------------------------------------------------------------------

if( it <= n_time_steps )then
    write(6,'(//8x, A, 1x,I2,3(1x,E15.7) )') &
          'prt4: icff, left, right, tree_eval ', &
           icff, left_node_value, right_node_value, &
           tree_evaluation(i_function,i_tree)
    write(6,'(8x, A, A )') 'prt4: left_node_value_string  ', &
                            trim( left_node_value_string )
    write(6,'(8x, A, A )') 'prt4: right_node_value_string ', &
                            trim( right_node_value_string )

    write(6,'(/8x, A, 1x,I6,1x,A,5x, A )') 'prt4: i_tree, tree_evaluation_string', &
                      i_tree, ':', trim( tree_evaluation_string(i_function,i_tree) )

endif !  it <= n_time_steps

return

END subroutine print4
