subroutine print_trees( itype ) 

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


logical bad_node


integer,intent(in) :: itype 
integer :: i
integer :: nop
integer :: message_len

logical :: Lplot


!real(kind=8), allocatable, dimension(:) :: answer
!real(kind=8), allocatable, dimension(:) :: output_array

real (kind=8) :: ssum, ssum2, totobs, dff

real (kind=8) :: xcount

!----------------------------------------------------------------------------------------


! print trees 



if( itype == 0 ) then
    write(6,'(A/)')  '0: trees before call to GP_Clean_Tree_Nodes'
else
    write(6,'(A/)')  '0: trees after call to GP_Clean_Tree_Nodes'
endif ! itype 

do  i_GP_individual = 1, n_GP_individuals 

    write(6,'(/A,1x,I6)')  '0: i_GP_indiv ', i_GP_individual
    write(6,'(A)') '0: i_tree nodes '
    write(6,'(A)') '            1  2  3  4  5  6  7  8  9 10 11 12 13 14 15'

    do  i_Tree=1,n_Trees

        write(6,'(I6,4x,20(1x,I2))' ) &
             i_tree, GP_Adult_Population_Node_Type(i_GP_individual, 1:n_nodes, i_tree) 

    enddo ! i_tree
enddo  ! i_GP_individual




return

end subroutine print_trees
