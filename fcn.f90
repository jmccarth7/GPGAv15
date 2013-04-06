subroutine fcn(mm,nn,x,fvec,iflag)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


use mpi
use mpi_module

use GP_parameters_module
use GP_variables_module
use GP_data_module

implicit none

integer(kind=4),intent(in)  :: mm  ! n_tsteps
integer(kind=4),intent(in)  :: nn  ! n_parms


real (kind=8) :: fvec(n_time_steps)

real (kind=8) :: x(n_parameters)
real (kind=8) :: sse


integer (kind=4) :: iflag


!real (kind=8) :: Data_Variance(n_CODE_equations)
!real (kind=8) :: ssum, ssum2, totobs, cff

!---------------------------------------------------------------------

! move the values you are trying to fit into 
! the initial conditions and variable terms


! set up the initial conditions

do i_CODE_equation=1,n_CODE_equations

  Runge_Kutta_Solution(0,i_CODE_equation)=dabs(x(i_CODE_equation))

  if( isnan( Runge_Kutta_Solution(0,i_CODE_equation) ) .or. &
        abs( Runge_Kutta_Solution(0,i_CODE_equation) )  > 1.0D20  )then

      write(6,'(A,1x,I6)') &
            'fcn: bad initial condition i_CODE_equation = ', &
                                        i_CODE_equation

      L_bad_result = .TRUE.
      iflag = -1
      return

  endif  ! isnan

enddo !  i_CODE_equation


! set the node_parameters array from the parameter array


Runge_Kutta_Node_Parameters = 0.

i_parameter=n_CODE_equations

do i_tree=1,n_trees
  do i_node=1,n_nodes

    if( GP_Individual_Node_Type(i_node,i_tree) .eq. 0) then  ! set the node_parameter

      i_parameter=i_parameter+1
      Runge_Kutta_Node_Parameters(i_node,i_tree)=dabs(x(i_parameter))

      if( isnan( Runge_Kutta_Node_Parameters(i_node,i_tree) )  .or. &
            abs( Runge_Kutta_Node_Parameters(i_node,i_tree) ) > 1.0D20 ) then

          write(6,'(A,2(1x,I6))') &
                'fcn: bad  Runge_Kutta_Node_Parameters for i_node, i_tree ', &
                                                           i_node, i_tree
          L_bad_result = .TRUE.
          iflag = -1
          return

      endif  ! isnan

      write(6,'(A,3(1x,I6),1x,E15.7)') &
            'fcn: i_node, i_tree, i_parameter, Runge_Kutta_Node_Parameters(i_node,i_tree) ', &
                  i_node, i_tree, i_parameter, Runge_Kutta_Node_Parameters(i_node,i_tree) 

    endif !  GP_individual_node_type(i_node,i_tree) .eq. 0

  enddo ! i_node
enddo  ! i_tree


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!  Runge_Kutta_Box_Model runs the RK process using the parameters
!  set above


L_bad_result = .FALSE.

call Runge_Kutta_Box_Model


if( L_bad_result ) then
    write(6,'(A,1x,I6,4x,L1)') &
          'fcn: aft call Runge_Kutta_Box_Model  myid, L_bad_result = ', &
                                                myid, L_bad_result
    iflag = -1
    return
endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! if the result of the RK process was good, compute the fvec (and maybe SSE)

!SSE=0.0D+0
do i_time_step=1,n_time_steps

  fvec(i_time_step)=0.

  !write(6,'(/A,1x,I6, 1x,I10)')'fcn: myid, i_time_step ', myid, i_time_step

  do i_CODE_equation=1,n_CODE_equations

    !if( abs(  data_variance(i_CODE_equation) ) > 1.0D-30  )then

        !write(6,'(A,2(1x,I6), 3(1x,E15.7))') &
        !      'fcn: myid, i_eqn, RK_soln, data_array, var ', &
        !            myid, i_CODE_equation,                   &
        !            Runge_Kutta_Solution(i_time_step,i_CODE_equation), &
        !            Data_Array(i_time_step,i_CODE_equation), &
        !            data_variance(i_CODE_equation)

        !write(6,'(A,2(1x,I6), 1x,E15.7)') &
        !      'fcn: myid, i_eqn, data_variance ', &
        !            myid, i_CODE_equation, data_variance(i_CODE_equation)

        fvec(i_time_step) = fvec(i_time_step)  +                                  &
            (   Data_Array(i_time_step,i_CODE_equation) -                         &
                  Runge_Kutta_Solution(i_time_step,i_CODE_equation)   )**2  /     &
                                                        Data_Variance(i_CODE_equation)

    !else

        !write(6,'(/A,1x,I6, 1x,I10)')'fcn: bad variance myid, i_time_step ', &
        !                                                myid, i_time_step
        !write(6,'(A,2(1x,I6), 1x,E15.7)') &
        !      'fcn: myid, i_CODE_equation, data_variance(i_CODE_equation)              ', &
        !            myid, i_CODE_equation, data_variance(i_CODE_equation)

    !endif ! abs(  data_variance(i_CODE_equation) ) > 1.0D-30



  enddo ! i_CODE_equation


  !SSE = SSE + fvec(i_time_step)

  !write(6,'(A,1x,I6, 1x,I6, 1x, E15.7)')&
  !      'fcn: myid, i_time_step, fvec ', &
  !            myid, i_time_step, fvec(i_time_step)

enddo ! i_time_step

!write(6,'(A,1x,I6,2x,E24.16)') 'fcn: myid, SSE = ',myid, SSE


return
end subroutine fcn
!234567890123456789012345678901234567890123456789012345678901234567890
