subroutine print_trees( n_indiv_start, n_indiv_stop, &
                        tree_type, tree_descrip )

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


integer,intent(in) :: n_indiv_start
integer,intent(in) :: n_indiv_stop
character(*),intent(in) :: tree_descrip

character(2),dimension( 1:n_nodes ) :: tree_type_string


integer(kind=4), intent(in), &
        dimension( 1:n_nodes, 1:n_trees, n_GP_individuals) :: tree_type

integer(kind=4) :: i_GP_individual
integer(kind=4) :: j_GP_individual
integer(kind=4) :: i_Tree
integer(kind=4) :: i_node
integer(kind=4) :: jj

!----------------------------------------------------------------------------------------


! print trees

write(GP_print_unit,'(/A)')  &
      'pt: ############################################################################'
write(GP_print_unit,'(A)')  tree_descrip

do  i_GP_individual = n_indiv_start, n_indiv_stop

    write(GP_print_unit,'(/A,1x,I6)')  'pt: i_GP_indiv ', i_GP_individual
    write(GP_print_unit,'(A)') 'pt: i_tree                    nodes '
    write(GP_print_unit,'(A)') '            1  2  3  4  5  6  7  8  9 10 11 12 13 14 15'

    j_GP_individual = i_GP_individual
    if( n_indiv_start ==  n_indiv_stop ) j_GP_individual = 1

    do  i_Tree=1,n_Trees

        !write(GP_print_unit,'(I6,4x,20(1x,I2))' ) &
        !     i_tree, Tree_Type( 1:n_nodes, i_tree,j_GP_individual)
        !write(GP_print_unit,'(I6,4x,20(1x,I2))' ) &

        tree_type_string = '  '
        do  i_node = 1, n_nodes

            !write(GP_print_unit,'(A,3(1x,I6))' ) &
            !     'pt: i_tree, i_node, Tree_Type( i_node, i_tree,j_GP_individual ) ', &
            !          i_tree, i_node, Tree_Type( i_node, i_tree,j_GP_individual )

            if( Tree_Type( i_node, i_tree, j_GP_individual) == -9999 )then
                tree_type_string(i_node)  = '**'
            else
                write(tree_type_string(i_node), '(I2)') &
                      tree_type( i_node, i_tree, j_GP_individual )
            endif ! tree_type(...-9999

            !write(GP_print_unit,'(A,2(1x,I6), 1x, A)' ) &
            !     'pt: i_tree, i_node, Tree_Type_string( i_node ) ', &
            !          i_tree, i_node, Tree_Type_string( i_node )

        enddo ! i_node


        !write(GP_print_unit,'(I6,4x,20(1x,I5))' ) &
        !     i_tree, ( Tree_Type(jj, i_tree,j_GP_individual), jj = 1, n_nodes )
        write(GP_print_unit,'(I6,4x,20(1x,A))' ) &
             i_tree, Tree_Type_string(1:n_nodes)

    enddo ! i_tree
enddo  ! i_GP_individual

write(GP_print_unit,'(A/)')  &
      'pt: ############################################################################'


return

end subroutine print_trees
