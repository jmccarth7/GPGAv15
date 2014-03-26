subroutine set_answer_arrays( buildTrees )

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
use GP_variables_module
use GA_Variables_module
use GP_Data_module


use Tree_Node_Factory_module
use class_Tree_Node


implicit none


integer(kind=4) :: i_GP_individual
integer(kind=4) :: i_tree
integer(kind=4) :: i_node
integer(kind=4) :: ii

logical :: buildTrees


!------------------------------------------------------------------------------


!Numerical_CODE_Solution(0,1:n_CODE_equations) = Runge_Kutta_Initial_Conditions ! Array Assignment
!Numerical_CODE_Initial_Conditions             = Runge_Kutta_Initial_Conditions ! Array Assignment


!!!allocate( GP_Trees(n_trees,1)  )

write(6,'(/A/)') 'saa: call Initialize_Model  '

call Initialize_Model( .true., .true., 6 )

write(6,'(/A/)') 'saa: aft call Initialize_Model  '


!------------------------------------------------------------------------------

!off include 'Lotka_Volterra_Example_Set_Up.f901'  ! replaced with routine init_values
!off include 'Franks_etal_NPZ_Mar_Bio_Example_Set_Up.f901'

! set the desired 'twin experiment' population node type
! and parameter using the info from the set up file

GP_Node_Type_Answer       = GP_Individual_Node_Type       ! Matrix Operation
GP_Node_Parameters_Answer = GP_Individual_Node_Parameters ! Matrix Operation


!--------------------------------------------------------------------------

if( L_unit50_output )then

    ! calculate array for writing on unit50.txt ( unit_gp_out )

    do  i_GP_Individual=1,n_GP_individuals
        GP_Node_Type_for_Plotting(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
                              GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)
    enddo

    if( myid == 0 )then
        write(unit_gp_out) GP_Node_Type_for_Plotting
    endif ! myid == 0

endif ! L_unit50_output

!--------------------------------------------------------------------------

! set the initial population node type using the info obtained
! from the set up file
! set the Initial Conditions, Model Parameters and Node Type
! for the 'twin experiment case'



! initialize the biological data fields

!Runge_Kutta_Solution(0,1:n_CODE_equations)    = Runge_Kutta_Initial_Conditions    ! Array Assignment


!Numerical_CODE_Solution(0,1:n_CODE_equations) = Runge_Kutta_Initial_Conditions ! Array Assignment
!Numerical_CODE_Initial_Conditions             = Runge_Kutta_Initial_Conditions ! Array Assignment

Numerical_CODE_Solution(0,1:n_CODE_equations) = Numerical_CODE_Initial_Conditions ! Array Assignment

!Runge_Kutta_Node_Parameters = GP_Individual_Node_Parameters  ! Matrix Operation
!Runge_Kutta_Node_Type       = GP_Individual_Node_Type        ! Matrix Operation


if( myid == 0 )then

    write(6,'(A)') ' '

    do  ii = 1, n_CODE_equations
        write(6,'(A,1x,I6,1x,E15.7)') 'saa: ii, Numerical_CODE_Initial_Conditions(ii) ', &
                                            ii, Numerical_CODE_Initial_Conditions(ii)
    enddo ! ii

    write(6,'(A)') ' '

    do  ii = 1, n_CODE_equations
        write(6,'(A,1x,I6,1x,E15.7)') 'saa: ii, Numerical_CODE_Solution(0,ii) ', &
                                            ii, Numerical_CODE_Solution(0,ii)
    enddo ! ii


    write(6,'(A)') ' '
    write(6,'(A,2(1x,I6))') 'saa: n_trees, n_nodes ', n_trees, n_nodes

!    write(6,'(/A)') &
!          'saa: i_tree  i_node  Runge_Kutta_Node_Parameters( i_node, i_tree ) '
!    do  i_tree = 1, n_trees
!        do  i_node = 1, n_nodes
!
!            if( Runge_Kutta_Node_Type( i_node, i_tree ) == 0     )then
!                write(6,'(2(1x,I8),6x,E15.7)') &
!                      i_tree, i_node, Runge_Kutta_Node_Parameters( i_node, i_tree )
!            endif ! Runge_Kutta_Node_Type( i_node, i_tree ) == 0
!
!        enddo ! i_node
!    enddo ! i_tree

!    write(6,'(//A)') &
!          'saa: i_tree  i_node  Runge_Kutta_Node_Type( i_node, i_tree ) '
!
!    do  i_tree = 1, n_trees
!        do  i_node = 1, n_nodes
!
!            if( Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999 )then
!                write(6,'(3(1x,I8))') &
!                        i_tree, i_node, Runge_Kutta_Node_Type( i_node, i_tree )
!            endif ! Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999
!
!        enddo ! i_node
!    enddo ! i_tree
!
!    write(6,'(A)') ' '

endif ! myid == 0



!------------------------------------------------------------------------


! run the Runge-Kutta model only once with proc 0

if( myid == 0 )then

    ! Runge_Kutta_Box_Model now put the time series in Numerical_CODE_Solution

    call Runge_Kutta_Box_Model( .false. )  ! 1 )

    !!Runge_Kutta_Solution = Numerical_CODE_Solution


endif ! myid == 0


return

end subroutine set_answer_arrays
