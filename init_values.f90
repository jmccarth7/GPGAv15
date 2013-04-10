subroutine init_values( icall  )


use mpi
use mpi_module
use GP_parameters_module
use GP_variables_module
use Runge_Kutta_Variables_module

implicit none


integer,intent(in)  :: icall




!-------------------------------------------------------------------------

if( myid == 0 )then
    write(GP_print_unit,'(/A,1x,A)')  'iv: model ', trim(model)
    write(GP_print_unit,'(A,1x,I6/)') 'iv: icall ', icall 
endif ! myid == 0 


if( trim(model) == 'NPZ' ) then

    call init_values_NPZ( icall )  
    if( icall == 0 ) return

elseif( trim(model) == 'LV' )then

    call init_values_LV( icall ) 
    if( icall == 0 ) return

endif ! trim(model) == 'NPZ' 


!----------------------------------------------------------------------------------

! print node types

!if( myid == 0 ) then
!
!    do  i_tree = 1, n_trees
!
!        write(GP_print_unit,'(/A,1x,i2)') 'i_tree = ', i_tree
!        write(GP_print_unit,'(32(1x,i2))') GP_Individual_node_type(1:n_nodes, i_tree)
!
!    enddo
!
!    write(GP_print_unit,'(/A)') ' '
!
!
!    !-------------------------------------------------------------------------------
!
!
!
!    do  i_node = 1,n_nodes
!        do  i_tree = 1,n_trees
!
!            if( GP_Individual_Node_Parameters(i_node,i_tree) > 0.0d0 )then
!
!                write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
!                      'iv: i_node, i_tree, GP_Individual_Node_Parameters(i_node,i_tree)  ', &
!                           i_node, i_tree, GP_Individual_Node_Parameters(i_node,i_tree)
!
!            endif
!        enddo ! i_tree
!
!    enddo ! i_node
!
!endif !  myid == 0

return

END subroutine init_values
