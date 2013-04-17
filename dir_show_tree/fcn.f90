subroutine fcn(mm,nn,x,fvec,iflag)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


use GP_parameters_module
use GP_variables_module
use GP_data_module

implicit none

integer(kind=4),intent(in)  :: mm  ! n_tsteps
integer(kind=4),intent(in)  :: nn  ! n_parms


real (kind=8) :: fvec(n_time_steps)

real (kind=8) :: x(n_maximum_number_parameters)

!real (kind=8) :: sse


integer (kind=4) :: iflag



!---------------------------------------------------------------------

! move the values you are trying to fit into
! the initial conditions and variable terms


! set up the initial conditions

do i_CODE_equation=1,n_CODE_equations

  Runge_Kutta_Solution(0,i_CODE_equation)=dabs(x(i_CODE_equation))

enddo !  i_CODE_equation


! set the node_parameters array from the parameter array


Runge_Kutta_Node_Parameters = 0.

i_parameter=n_CODE_equations

do i_tree=1,n_trees
  do i_node=1,n_nodes

    if( GP_Individual_Node_Type(i_node,i_tree) .eq. 0) then  ! set the node_parameter

      i_parameter=i_parameter+1
      Runge_Kutta_Node_Parameters(i_node,i_tree)=dabs(x(i_parameter))


      !write(6,'(A,3(1x,I6),1x,E15.7)') &
      !      'fcn: i_node, i_tree, i_parameter, Runge_Kutta_Node_Parameters(i_node,i_tree) ', &
      !            i_node, i_tree, i_parameter, Runge_Kutta_Node_Parameters(i_node,i_tree)

    endif !  GP_individual_node_type(i_node,i_tree) .eq. 0

  enddo ! i_node
enddo  ! i_tree


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!  Runge_Kutta_Box_Model runs the RK process using the parameters
!  set above


L_bad_result = .FALSE.

call Runge_Kutta_Box_Model

!write(6,'(A)') 'fcn: after call Runge_Kutta_Box_Model '                                                



return
end subroutine fcn
!234567890123456789012345678901234567890123456789012345678901234567890
