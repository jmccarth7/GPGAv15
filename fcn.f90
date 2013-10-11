subroutine fcn(mm,nn,x,fvec,iflag)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


use mpi
use mpi_module

use GA_parameters_module
use GP_parameters_module
use GP_variables_module
use GP_data_module

implicit none

integer(kind=4),intent(in)  :: mm  ! n_tsteps
integer(kind=4),intent(in)  :: nn  ! n_parms


real(kind=8) :: fvec(n_time_steps)

!real(kind=8) :: x(n_maximum_number_parameters)
real(kind=8) :: x( nn )

real(kind=8) :: sse_local 

integer(kind=4) :: iflag

integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node

integer(kind=4) :: i_CODE_equation
integer(kind=4) :: i_time_step
integer(kind=4) :: i_parameter

logical,parameter :: L_GP_print = .TRUE.

!---------------------------------------------------------------------

!if( GP_para_flag .and. myid == 1 )then
!    write(GP_print_unit,'(A,2(1x,I6))') &
!           'fcn: n_CODE_equations, nn', &
!                 n_CODE_equations, nn
!endif ! myid == 1

! move the values you are trying to fit into
! the initial conditions and variable terms


! set up the initial conditions

Runge_Kutta_Solution = 0.0D0

do i_CODE_equation=1,n_CODE_equations

  Runge_Kutta_Solution(0,i_CODE_equation)=dabs(x(i_CODE_equation))

  if( isnan( Runge_Kutta_Solution(0,i_CODE_equation) ) .or. &
        abs( Runge_Kutta_Solution(0,i_CODE_equation) )  > 1.0D20  )then

      !if( L_GP_print )then
      !    write(GP_print_unit,'(A,1x,I6)') &
      !          'fcn: bad initial condition i_CODE_equation = ', &
      !                                      i_CODE_equation
      !endif ! L_GP_print

      L_bad_result = .TRUE.
      iflag = -1
      return

  endif  ! isnan

  !if( GP_para_flag .and. myid == 1 )then
  !    write(GP_print_unit,'(A,1x,I6, 2(1x,E15.7))') &
  !     'fcn: i_CODE_equation, Runge_Kutta_Solution(0,), X', &
  !           i_CODE_equation, Runge_Kutta_Solution(0,i_CODE_equation),&
  !                                                 x(i_CODE_equation)
  !endif ! myid == 1
enddo !  i_CODE_equation


! set the node_parameters array from the parameter array


Runge_Kutta_Node_Parameters = 0.0D0
Runge_Kutta_Node_Type       = -9999    ! jjm 20130417

i_parameter = n_CODE_equations

do i_tree=1,n_trees
  do i_node=1,n_nodes

    Runge_Kutta_Node_Type(i_node,i_tree) = &
            GP_Individual_Node_Type(i_node,i_tree) ! jjm 20130417

    !if( GP_para_flag .and. myid == 1 )then
    !    if( GP_Individual_Node_Type(i_node,i_tree) > -9999 )then
    !        write(GP_print_unit,'(A,4(1x,I6))') &
    !           'fcn: myid, i_tree, i_node, GP_Individual_Node_Type', &
    !                 myid, i_tree, i_node, GP_Individual_Node_Type(i_node,i_tree)
    !    endif ! GP_Individual_Node_Type(i_node,i_tree) > -9999 
    !endif !  myid == 1

    if( GP_Individual_Node_Type(i_node,i_tree) .eq. 0) then  ! set the node_parameter

      i_parameter=i_parameter+1
      Runge_Kutta_Node_Parameters(i_node,i_tree)=dabs(x(i_parameter))

      if( isnan( Runge_Kutta_Node_Parameters(i_node,i_tree) )  .or. &
            abs( Runge_Kutta_Node_Parameters(i_node,i_tree) ) > 1.0D20 ) then

          !if( L_GP_print )then
          !    write(GP_print_unit,'(A,2(1x,I6))') &
          !          'fcn: bad  Runge_Kutta_Node_Parameters for i_tree, i_node ', &
          !                                                     i_tree, i_node
          !endif ! L_GP_print

          L_bad_result = .TRUE.
          iflag = -1
          return

      endif  ! isnan

      !if( L_GP_print )then
      !    if( GP_para_flag .and. myid == 1 )then
      !        write(GP_print_unit,'(A,5(1x,I4),1x,E15.7)') &
      !        'fcn: myid, i_tree, i_node, i_parameter, nn,  Runge_Kutta_Node_Params', &
      !              myid, i_tree, i_node, i_parameter, nn,  &
      !                                Runge_Kutta_Node_Parameters(i_node,i_tree)
      !    endif !  myid == 1
      !endif ! L_GP_print

    endif !  GP_individual_node_type(i_node,i_tree) .eq. 0

  enddo ! i_node
enddo  ! i_tree

!if( L_GP_print )then
!    if( GP_para_flag .and. myid == 1 )then
!        write(6,'(/A,3(1x,I4)/)') &
!               'fcn: myid, i_parameter, nn', myid, i_parameter, nn
!    endif !  myid == 1
!endif ! L_GP_print

!---------------------------------------------------------------------------------

! jjm 20130417 >>>>>>>>>>>>>>>
!if( L_GP_print )then
!    do i_tree=1,n_trees
!      do i_node=1,n_nodes
!          if( GP_para_flag .and. myid == 1) then
!              if( Runge_Kutta_Node_Type( i_node, i_tree ) == 0   )then
!                  write(GP_print_unit,'(A,3(1x,I6),1x,E15.7)') &
!                   'fcn: myid, i_node, i_tree, Runge_Kutta_Node_Parameters(i_node,i_tree) ', &
!                         myid, i_node, i_tree, Runge_Kutta_Node_Parameters(i_node,i_tree)
!              endif ! Runge_Kutta_Node_Type( i_node, i_tree ) == 0
!          endif ! myid == 1   NOTE:  1
!      enddo ! i_node
!    enddo  ! i_tree
!    do i_tree=1,n_trees
!      do i_node=1,n_nodes
!          if( GP_para_flag .and. myid == 1) then
!              if( Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999  )then
!                  write(GP_print_unit,'(A,4(1x,I6))') &
!                   'fcn: myid, i_node, i_tree, Runge_Kutta_Node_Type(i_node,i_tree) ', &
!                         myid, i_node, i_tree, Runge_Kutta_Node_Type(i_node,i_tree)
!              endif ! Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999
!          endif ! myid == 1   NOTE:  1
!      enddo ! i_node
!    enddo  ! i_tree
!endif ! L_GP_print
! jjm 20130417 <<<<<<<<<<<<<<<

!---------------------------------------------------------------------------------


!  Runge_Kutta_Box_Model runs the RK process using the parameters
!  set above


L_bad_result = .FALSE.

!if( L_GP_print .and. GP_para_flag .and. myid == 1 )then
!    write(GP_print_unit,'(A,1x,I6)') &
!          'fcn: call Runge_Kutta_Box_Model  myid', &
!                                            myid
!endif ! L_GP_print


call Runge_Kutta_Box_Model


if( L_bad_result ) then
    !if( L_GP_print ) then !.and. GP_para_flag .and. myid == 1 )then
    !    write(GP_print_unit,'(A,1x,I6,4x,L1)') &
    !          'fcn: aft call Runge_Kutta_Box_Model  myid, L_bad_result = ', &
    !                                                myid, L_bad_result
    !endif ! L_GP_print
    iflag = -1
    return
endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! if the result of the RK process was good, compute the fvec (and maybe sse_local)

!sse_local=0.0D+0
do i_time_step=1,n_time_steps

  fvec(i_time_step)=0.0D0

  !if( L_GP_print .and. GP_para_flag .and. myid == 1 .and. i_time_step == n_time_steps )then
  !    write(GP_print_unit,'(/A,1x,I6, 1x,I10)') &
  !      'fcn: myid, i_time_step', myid, i_time_step
  !endif ! L_GP_print

  do  i_CODE_equation=1,n_CODE_equations

      !if( L_GP_print .and. GP_para_flag .and. myid == 1 .and. i_time_step == n_time_steps )then
      !    write(GP_print_unit,'(A,2(1x,I6), 3(1x,E15.7))') &
      !          'fcn: myid, i_eqn, RK_soln, data_array, var_inv ', &
      !                myid, i_CODE_equation,                   &
      !                Runge_Kutta_Solution(i_time_step,i_CODE_equation), &
      !                Data_Array(i_time_step,i_CODE_equation), &
      !                data_variance_inv(i_CODE_equation)
      !    write(GP_print_unit,'(A,2(1x,I6), 1x,E15.7)') &
      !          'fcn: myid, i_eqn, data_variance_inv ', &
      !                myid, i_CODE_equation, data_variance_inv(i_CODE_equation)
      !endif ! L_GP_print


      !fvec(i_time_step) = fvec(i_time_step)  +                             &
      ! (   Data_Array(i_time_step,i_CODE_equation) -                       &
      !     Runge_Kutta_Solution(i_time_step,i_CODE_equation)   )**2  /     &
      !                                     Data_Variance(i_CODE_equation)

      fvec(i_time_step) = fvec(i_time_step)  +                                &
          (   Data_Array(i_time_step,i_CODE_equation) -                       &
              Runge_Kutta_Solution(i_time_step,i_CODE_equation)   )**2  *     &
                                            Data_Variance_inv(i_CODE_equation)

  enddo ! i_CODE_equation


  !sse_local = sse_local + fvec(i_time_step)

  !if( L_GP_print .and. GP_para_flag .and. myid == 1 .and. i_time_step == n_time_steps )then
  !   write(GP_print_unit,'(A,1x,I6, 1x,I6, 1x, E15.7)')&
  !          'fcn: myid, i_time_step, fvec', &
  !                myid, i_time_step, fvec(i_time_step)
  !endif ! L_GP_print

enddo ! i_time_step

!if( L_GP_print .and. GP_para_flag .and. myid == 1  )then
!    write(GP_print_unit,'(A,1x,I6,2x,E24.16)') 'fcn: myid, sse_local = ',myid, sse_local
!endif ! L_GP_print


return


end subroutine fcn
