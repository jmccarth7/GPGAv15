subroutine print_trees( n_indiv_start, n_indiv_stop, tree_type, tree_descrip )

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

!integer :: i

integer(kind=4), intent(in), &
        dimension(n_GP_individuals, 1:n_nodes, 1:n_trees) :: tree_type

!----------------------------------------------------------------------------------------


! print trees

write(GP_print_unit,'(/A/)')  &
      'pt: ############################################################################'
write(GP_print_unit,'(A/)')  tree_descrip

do  i_GP_individual = n_indiv_start, n_indiv_stop

    write(GP_print_unit,'(/A,1x,I6)')  'pt: i_GP_indiv ', i_GP_individual
    write(GP_print_unit,'(A)') 'pt: i_tree                    nodes '
    write(GP_print_unit,'(A)') '            1  2  3  4  5  6  7  8  9 10 11 12 13 14 15'

    do  i_Tree=1,n_Trees

        write(GP_print_unit,'(I6,4x,20(1x,I2))' ) &
             i_tree, Tree_Type(i_GP_individual, 1:n_nodes, i_tree)

    enddo ! i_tree
enddo  ! i_GP_individual

write(GP_print_unit,'(/A/)')  &
      'pt: ############################################################################'


return

end subroutine print_trees
