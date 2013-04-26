subroutine set_answer_arrays( )

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


integer :: i
integer :: nop
integer :: message_len

logical :: Lplot


!real(kind=8), allocatable, dimension(:) :: answer
!real(kind=8), allocatable, dimension(:) :: output_array

real (kind=8) :: ssum, ssum2, totobs, dff

real (kind=8) :: xcount

!------------------------------------------------------------------------------



!------------------------------------------------------------------

!off include 'Lotka_Volterra_Example_Set_Up.f901'  ! replaced with routine init_values
!off include 'Franks_etal_NPZ_Mar_Bio_Example_Set_Up.f901'

! set the desired 'twin experiment' population node type 
! and parameter using the info from the set up file

GP_Node_Type_Answer       = GP_Individual_Node_Type       ! Matrix Operation
GP_Node_Parameters_Answer = GP_Individual_Node_Parameters ! Matrix Operation


do i_GP_Individual=1,n_GP_individuals
  GP_Node_Type_for_Plotting(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
                        GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)
enddo

if( myid == 0 )then
    write(unit_gp_out) GP_Node_Type_for_Plotting
endif ! myid == 0

! set the initial population node type using the info obtained 
! from the set up file
! set the Initial Conditions, Model Parameters and Node Type 
! for the 'twin experiment case'

GP_Population_Node_Type(1,1:n_nodes,1:n_trees)=GP_Individual_Node_Type(1:n_nodes,1:n_trees)


! initialize the biological data fields

Runge_Kutta_Solution(0,1:n_CODE_equations)=Runge_Kutta_Initial_Conditions ! Array Assignment

Runge_Kutta_Node_Parameters = GP_Individual_Node_Parameters  ! Matrix Operation
Runge_Kutta_Node_Type=GP_Individual_Node_Type                ! Matrix Operation


!if( myid == 0 )then
!    !write(6,'(A)') ' '
!
!    do  i_tree = 1, n_trees
!        do  i_node = 1, n_nodes
!
!            if( Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999 )then 
!                !write(6,'(A,2(1x,I6),1x,E15.7)') &
!                !  'saa: i_tree, i_node, Runge_Kutta_Node_Parameters( i_node, i_tree ) ', &
!                !        i_tree, i_node, Runge_Kutta_Node_Parameters( i_node, i_tree )
!            endif ! Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999 
!
!        enddo ! i_node 
!    enddo ! i_tree
!
!    !write(6,'(A)') ' '
!
!    do  i_tree = 1, n_trees
!        do  i_node = 1, n_nodes
!
!            if( Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999 )then 
!                !write(6,'(A,3(1x,I6))') &
!                !  'saa: i_tree, i_node, Runge_Kutta_Node_Type( i_node, i_tree ) ', &
!                !        i_tree, i_node, Runge_Kutta_Node_Type( i_node, i_tree )
!            endif ! Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999 
!
!        enddo ! i_node 
!    enddo ! i_tree
!
!    write(6,'(A)') ' '
!
!endif ! myid == 0



!------------------------------------------------------------------------


! run the Runge-Kutta model only once with proc 0

if( myid == 0 )then
    call Runge_Kutta_Box_Model
endif ! myid == 0


return

end subroutine set_answer_arrays
