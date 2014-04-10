subroutine print_debug_integer_node_tree( iunit, label, input_array  )




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


integer(kind=4),intent(in) :: iunit

integer(kind=4) :: i_GP_Individual

integer(kind=4) :: i_tree
integer(kind=4) :: i_node


integer(kind=4) :: ierr


integer(kind=4), dimension( 1:n_Nodes,1:n_Trees, 1:n_GP_individuals ) :: &
                         input_array

!--------------------------------------------------------------------------------

write(iunit,'(/A)') 'pd3: entry print_debug3'


!! debug
write(iunit,'(/A,1x,A)') 'pd3: print ', trim(label)
write(iunit,'(/A)') &
   'pd3: i_GP_Individual  i_tree  i_node   input_array(i_Node,i_Tree, i_GP_individual )'

do  i_GP_individual = 1, n_GP_individuals
    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes

            if( input_array(i_Node,i_Tree, i_GP_individual ) > -9999 )then

                write(iunit,'(3x,4(1x,I10))',iostat = ierr ) &
                  i_GP_Individual, i_tree, i_node, &
                         input_array(i_Node,i_Tree, i_GP_individual )

                if( ierr /= 0 )then
                    write(iunit,*) 'pd3: write error  ierr = ', ierr
                endif ! ierr /= 0

            endif ! input_array(i_Node,i_Tree, i_GP_individual ) > 0.0d0

        enddo
    enddo
enddo ! i_GP_individual



!-------------------------------------------------------------------------------------------------

write(iunit,'(/A/)') 'pd3: at return   '

return

end subroutine print_debug_integer_node_tree
