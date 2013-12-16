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

real(kind=8) :: x(n_maximum_number_parameters)
!real(kind=8) :: x( nn )

real(kind=8) :: sse_local

integer(kind=4) :: iflag

integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node
integer(kind=4) :: ii

integer(kind=4) :: i_CODE_equation
integer(kind=4) :: i_time_step
integer(kind=4) :: i_parameter

logical,parameter :: L_GP_print = .TRUE.

!---------------------------------------------------------------------


!write(6,'(A,1x,I5,4x,L1)') 'fcn: myid, L_GA_print   ', myid, L_GA_print 
!write(6,'(A,1x,I5,1x,I5)') 'fcn: myid, GA_print_unit', myid, GA_print_unit

!if( GP_para_flag .and. myid == 3 )then
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

  Numerical_CODE_Solution(0,i_CODE_equation) = dabs( x(i_CODE_equation) )

  !if( GP_para_flag )then
  !    write(GP_print_unit,'(A,1x,I6, 2(1x,E15.7))') &
  !     'fcn: i_CODE_equation, Runge_Kutta_Solution(0,), X', &
  !           i_CODE_equation, Runge_Kutta_Solution(0,i_CODE_equation),&
  !                                                 x(i_CODE_equation)
  !endif ! GP_para_flag

  !if( L_GA_print )then
  !    write(6,'(A,1x,I6, 2(1x,E15.7))') &
  !     'fcn: i_CODE_equation, Runge_Kutta_Solution(0,), X', &
  !           i_CODE_equation, Runge_Kutta_Solution(0,i_CODE_equation),&
  !                                                 x(i_CODE_equation)
  !    write(GA_print_unit,'(A,1x,I6, 2(1x,E15.7))') &
  !     'fcn: i_CODE_equation, Runge_Kutta_Solution(0,), X', &
  !           i_CODE_equation, Runge_Kutta_Solution(0,i_CODE_equation),&
  !                                                 x(i_CODE_equation)
  !    write(6,'(A,1x,I6, 2(1x,E15.7))') &
  !     'fcn: i_CODE_equation, Numerical_CODE_Solution(0,), X', &
  !           i_CODE_equation, Numerical_CODE_Solution(0,i_CODE_equation),&
  !                                                    x(i_CODE_equation)
  !    write(GA_print_unit,'(A,1x,I6, 2(1x,E15.7))') &
  !     'fcn: i_CODE_equation, Numerical_CODE_Solution(0,), X', &
  !           i_CODE_equation, Numerical_CODE_Solution(0,i_CODE_equation),&
  !                                                    x(i_CODE_equation)
  !endif ! L_GA_print

  if( isnan( Runge_Kutta_Solution(0,i_CODE_equation) ) .or. &
        abs( Runge_Kutta_Solution(0,i_CODE_equation) )  > 1.0D20  )then

      !if( L_GP_print )then
      !    write(GP_print_unit,'(A,2(1x,I6))') &
      !          'fcn: bad initial condition myid,i_CODE_equation = ', &
      !                                      myid,i_CODE_equation
      !endif ! L_GP_print

      L_bad_result = .TRUE.
      iflag = -1
      return

  endif  ! isnan


  if( isnan( Numerical_CODE_Solution(0,i_CODE_equation) ) .or. &
        abs( Numerical_CODE_Solution(0,i_CODE_equation) )  > 1.0D20  )then

      !if( L_GP_print )then
      !    write(GP_print_unit,'(A,2(1x,I6))') &
      !          'fcn: bad initial condition myid,i_CODE_equation = ', &
      !                                      myid,i_CODE_equation
      !endif ! L_GP_print

      L_bad_result = .TRUE.
      iflag = -1
      return

  endif  ! isnan

enddo !  i_CODE_equation


! set the node_parameters array from the parameter array

!-----------------------------------------------------------------

!Runge_Kutta_Node_Parameters = 0.0D0                        ! 20131209
!Runge_Kutta_Node_Type       = -9999    ! jjm 20130417

do  i_tree=1,n_trees                        ! 20131209
    do  i_node=1,n_nodes
        Runge_Kutta_Node_Parameters(i_node,i_tree) = 0.0D0
        Runge_Kutta_Node_Type(i_node,i_tree)       = -9999    ! jjm 20130417
  enddo ! i_node
enddo  ! i_tree

!-----------------------------------------------------------------


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
    ! if( L_GA_print )then
    !     if( GP_Individual_Node_Type(i_node,i_tree) > -9999 )then
    !         write(GA_print_unit,'(A,4(1x,I6))') &
    !            'fcn: myid, i_tree, i_node, GP_Individual_Node_Type', &
    !                  myid, i_tree, i_node, GP_Individual_Node_Type(i_node,i_tree)
    !     endif ! GP_Individual_Node_Type(i_node,i_tree) > -9999
    ! endif ! L_GA_print

    if( GP_Individual_Node_Type(i_node,i_tree) .eq. 0) then  ! set the node_parameter

      i_parameter=i_parameter+1
      Runge_Kutta_Node_Parameters(i_node,i_tree)   = dabs(x(i_parameter))
      GP_Individual_Node_Parameters(i_node,i_tree) = dabs(x(i_parameter))

      !if( L_GP_print )then
      !    if( GP_para_flag )then
      !        write(GP_print_unit,'(A,5(1x,I4),1x,E15.7)') &
      !        'fcn: myid, i_tree, i_node, i_parameter, nn,  Runge_Kutta_Node_Params', &
      !              myid, i_tree, i_node, i_parameter, nn,  &
      !                                Runge_Kutta_Node_Parameters(i_node,i_tree)
      !    endif !  GP_para_flag
      !endif ! L_GP_print
      !if( L_GA_print )then
      !    write(GA_print_unit,'(A,5(1x,I4),1x,E15.7)') &
      !    'fcn: myid, i_tree, i_node, i_parameter, nn,  Runge_Kutta_Node_Params', &
      !          myid, i_tree, i_node, i_parameter, nn,  &
      !                            Runge_Kutta_Node_Parameters(i_node,i_tree)
      !    write(GA_print_unit,'(A,5(1x,I4),1x,E15.7)') &
      !    'fcn: myid, i_tree, i_node, i_parameter, nn,  GP_Indiv_Node_Params', &
      !          myid, i_tree, i_node, i_parameter, nn,  &
      !                            GP_Individual_Node_Parameters(i_node,i_tree)
      !endif ! L_GA_print

      if( isnan( Runge_Kutta_Node_Parameters(i_node,i_tree) )  .or. &
            abs( Runge_Kutta_Node_Parameters(i_node,i_tree) ) > 1.0D20 ) then

          !if( L_ga_print )then
          !    write(GA_print_unit,'(A,3(1x,I6))') &
          !     'fcn: bad  Runge_Kutta_Node_Parameters for myid, i_tree, i_node ', &
          !                                                myid, i_tree, i_node
          !endif ! L_ga_print

          !if( L_GP_print )then
          !    write(GP_print_unit,'(A,3(1x,I6))') &
          !     'fcn: bad  Runge_Kutta_Node_Parameters for myid, i_tree, i_node ', &
          !                                                myid, i_tree, i_node
          !endif ! L_GP_print

          L_bad_result = .TRUE.
          iflag = -1
          return

      endif  ! isnan
      if( isnan( GP_Individual_Node_Parameters(i_node,i_tree) )  .or. &
            abs( GP_Individual_Node_Parameters(i_node,i_tree) ) > 1.0D20 ) then

          !if( L_ga_print )then
          !    write(GA_print_unit,'(A,3(1x,I6))') &
          !     'fcn: bad  GP_Individual_Node_Parameters for myid, i_tree, i_node ', &
          !                                                myid, i_tree, i_node
          !endif ! L_ga_print

          !if( L_GP_print )then
          !    write(GP_print_unit,'(A,3(1x,I6))') &
          !     'fcn: bad  GP_Individual_Node_Parameters for myid, i_tree, i_node ', &
          !                                                myid, i_tree, i_node
          !endif ! L_GP_print

          L_bad_result = .TRUE.
          iflag = -1
          return

      endif  ! isnan


      !if( L_ga_print )then
      !    write(GA_print_unit,'(A,4(1x,I4),1x,E15.7)') &
      !     'fcn: myid, i_tree, i_node, i_parameter, Runge_Kutta_Node_Params', &
      !           myid, i_tree, i_node, i_parameter, Runge_Kutta_Node_Parameters(i_node,i_tree)
      !    write(GA_print_unit,'(A,4(1x,I4),1x,E15.7)') &
      !     'fcn: myid, i_tree, i_node, i_parameter, GP_Indiv_Node_Params', &
      !           myid, i_tree, i_node, i_parameter, GP_Individual_Node_Parameters(i_node,i_tree)
      !endif ! L_ga_print

    endif !  GP_individual_node_type(i_node,i_tree) .eq. 0

  enddo ! i_node
enddo  ! i_tree

!if( L_GP_print )then
!    if( GP_para_flag .and. myid == 3 )then
!        write(6,'(/A,3(1x,I4)/)') &
!               'fcn: myid, i_parameter, nn', myid, i_parameter, nn
!    endif !  myid == 3
!endif ! L_GP_print

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,2(1x,I4)/)') &
!           'fcn: myid, i_parameter', myid, i_parameter
!endif ! L_ga_print

!---------------------------------------------------------------------------------

! jjm 20130417 >>>>>>>>>>>>>>>
!if( L_GP_print )then
!    do i_tree=1,n_trees
!      do i_node=1,n_nodes
!          if( GP_para_flag .and. myid == 1) then
!              if( Runge_Kutta_Node_Type( i_node, i_tree ) == 0   )then
!                  write(GP_print_unit,'(A,3(1x,I6),1x,E15.7)') &
!                   'fcn: myid, i_tree, i_node, Runge_Kutta_Node_Parameters(i_node,i_tree) ', &
!                         myid, i_tree, i_node, Runge_Kutta_Node_Parameters(i_node,i_tree)
!              endif ! Runge_Kutta_Node_Type( i_node, i_tree ) == 0
!          endif ! myid == 1   NOTE:  1
!      enddo ! i_node
!    enddo  ! i_tree
!    do i_tree=1,n_trees
!      do i_node=1,n_nodes
!          if( GP_para_flag .and. myid == 1) then
!              if( Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999  )then
!                  write(GP_print_unit,'(A,4(1x,I6))') &
!                   'fcn: myid, i_tree, i_node, Runge_Kutta_Node_Type(i_node,i_tree) ', &
!                         myid, i_tree, i_node, Runge_Kutta_Node_Type(i_node,i_tree)
!              endif ! Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999
!          endif ! myid == 1   NOTE:  1
!      enddo ! i_node
!    enddo  ! i_tree
!endif ! L_GP_print
!if( L_ga_print )then
!    do i_tree=1,n_trees
!      do i_node=1,n_nodes
!          if( myid == 1) then
!              if( Runge_Kutta_Node_Type( i_node, i_tree ) == 0   )then
!                  write(GA_print_unit,'(A,3(1x,I6),1x,E15.7)') &
!                  'fcn: myid, i_tree, i_node, Runge_Kutta_Node_Parameters(i_node,i_tree) ', &
!                        myid, i_tree, i_node, Runge_Kutta_Node_Parameters(i_node,i_tree)
!              endif ! Runge_Kutta_Node_Type( i_node, i_tree ) == 0
!          endif ! myid == 1   NOTE:  1
!      enddo ! i_node
!    enddo  ! i_tree
!    do i_tree=1,n_trees
!      do i_node=1,n_nodes
!          if( myid == 1) then
!              if( Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999  )then
!                  write(GA_print_unit,'(A,4(1x,I6))') &
!                   'fcn: myid, i_tree, i_node, Runge_Kutta_Node_Type(i_node,i_tree) ', &
!                         myid, i_tree, i_node, Runge_Kutta_Node_Type(i_node,i_tree)
!              endif ! Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999
!          endif ! myid == 1   NOTE:  1
!      enddo ! i_node
!    enddo  ! i_tree
!endif ! L_ga_print
! jjm 20130417 <<<<<<<<<<<<<<<

!---------------------------------------------------------------------------------

! set up the GP_Trees for the Runge_Kutta integration


! Initialize_Model calls build_trees which makes the GP_Trees

!if( L_ga_print )then
!    write(GA_print_unit,'(/A/)') 'fcn: call Initialize_Model(.true.)'
!endif ! L_ga_print 


! sets buildtrees = .true. in initialize_model

call Initialize_Model( .true. )   ! call build_trees

!call Initialize_Model(.false.)   ! call Deserialize_Trees

!if( L_ga_print )then
!    write(GA_print_unit,'(/A/)') 'fcn: aft call Initialize_Model(.true.)'
!endif ! L_ga_print 


!------------------------------------------------------------------------------

! Generate PDF representation of trees


!if( myid == 0 )then
    write(6,'(/A/)') 'fcn: call Generate_Dot_Graph'

    call Generate_Dot_Graph( GP_Trees(:,1), n_Trees, output_dir )

    !write(6,'(/A/)') 'fcn: aft call Generate_Dot_Graph'
!endif ! myid == 0


! Write trees to disk

!if( myid == 0 )then
    write(6,'(/A/)') 'fcn: call Serialize_Trees   '

    call Serialize_Trees( GP_Trees(:,:), n_Trees, n_Tracked_Resources, output_dir )

    !write(6,'(/A/)') 'fcn: aft call Serialize_Trees   '
!endif ! myid == 0

!------------------------------------------------------------------------------



! initialize the biological data fields


! Runge_Kutta_Solution set above from the x array
! Runge_Kutta_Solution(0,1:n_CODE_equations)    = Runge_Kutta_Initial_Conditions    ! Array Assignment


Numerical_CODE_Solution(0,1:n_CODE_equations) = Runge_Kutta_Solution(0,1:n_CODE_equations)
Numerical_CODE_Initial_Conditions             = Numerical_CODE_Solution(0,1:n_CODE_equations)


! set above
!Runge_Kutta_Node_Parameters = GP_Individual_Node_Parameters  ! Matrix Operation
!Runge_Kutta_Node_Type       = GP_Individual_Node_Type        ! Matrix Operation


!if( myid == 1 )then
!    write(GA_print_unit,'(A)') ' '
!
!    do  ii = 1, n_CODE_equations
!        write(GA_print_unit,'(A,1x,I6,1x,E15.7)') 'fcn: ii, Numerical_CODE_Initial_Conditions(ii)', &
!                                                        ii, Numerical_CODE_Initial_Conditions(ii)
!    enddo ! ii
!
!    write(GA_print_unit,'(A)') ' '
!
!    do  ii = 1, n_CODE_equations
!        write(GA_print_unit,'(A,1x,I6,1x,E15.7)') 'fcn: ii, Numerical_CODE_Solution(0,ii)        ', &
!                                                        ii, Numerical_CODE_Solution(0,ii)
!    enddo ! ii
!
!
!    write(GA_print_unit,'(A)') ' '
!
!    write(GA_print_unit,'(A,2(1x,I6))') 'fcn: n_trees, n_nodes ', n_trees, n_nodes
!
!    write(GA_print_unit,'(/A)') &
!          'fcn: i_tree  i_node  Runge_Kutta_Node_Parameters( i_node, i_tree ) '
!    do  i_tree = 1, n_trees
!        do  i_node = 1, n_nodes
!
!            if( Runge_Kutta_Node_Type( i_node, i_tree ) == 0     )then
!                write(GA_print_unit,'(2(1x,I8),6x,E15.7)') &
!                      i_tree, i_node, Runge_Kutta_Node_Parameters( i_node, i_tree )
!            endif ! Runge_Kutta_Node_Type( i_node, i_tree ) == 0
!
!        enddo ! i_node
!    enddo ! i_tree
!
!    write(GA_print_unit,'(/A)') &
!          'fcn: i_tree  i_node  GP_Individual_Node_Parameters( i_node, i_tree ) '
!    do  i_tree = 1, n_trees
!        do  i_node = 1, n_nodes
!
!            if( GP_Individual_Node_Type( i_node, i_tree ) == 0     )then
!                write(GA_print_unit,'(2(1x,I8),6x,E15.7)') &
!                      i_tree, i_node, GP_Individual_Node_Parameters( i_node, i_tree )
!            endif ! GP_Individual_Node_Type( i_node, i_tree ) == 0
!
!        enddo ! i_node
!    enddo ! i_tree
!
!    write(GA_print_unit,'(//A)') &
!          'fcn: i_tree  i_node  Runge_Kutta_Node_Type( i_node, i_tree ) '
!
!    do  i_tree = 1, n_trees
!        do  i_node = 1, n_nodes
!
!            if( Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999 )then
!                write(GA_print_unit,'(3(1x,I8))') &
!                        i_tree, i_node, Runge_Kutta_Node_Type( i_node, i_tree )
!            endif ! Runge_Kutta_Node_Type( i_node, i_tree ) /= -9999
!
!        enddo ! i_node
!    enddo ! i_tree
!
!    write(GA_print_unit,'(A)') ' '
!
!endif ! myid == 1



!---------------------------------------------------------------------------------

!  Runge_Kutta_Box_Model runs the RK process using the parameters
!  set above


L_bad_result = .FALSE.

!if( L_GP_print .and. GP_para_flag .and. myid == 1 )then
!    write(GP_print_unit,'(A,1x,I6)') &
!          'fcn: call Runge_Kutta_Box_Model  myid', &
!                                            myid
!endif ! L_GP_print
!if( L_ga_print )then
!    write(GA_print_unit,'(A,1x,I6)') &
!          'fcn: call Runge_Kutta_Box_Model  myid', &
!                                            myid
!endif ! L_ga_print


call Runge_Kutta_Box_Model


if( L_bad_result ) then

    !if( L_GP_print ) then !.and. GP_para_flag .and. myid == 1 )then
    !    write(GP_print_unit,'(A,1x,I6,4x,L1)') &
    !          'fcn: aft call Runge_Kutta_Box_Model  myid, L_bad_result = ', &
    !                                                myid, L_bad_result
    !endif ! L_GP_print
    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,1x,I6,4x,L1)') &
    !          'fcn: aft call Runge_Kutta_Box_Model  myid, L_bad_result = ', &
    !                                                myid, L_bad_result
    !endif ! L_ga_print
    iflag = -1
    return

endif ! L_bad_result


Runge_Kutta_Solution = Numerical_CODE_Solution


!stop ! debug only

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! if the result of the RK process was good, compute the fvec (and maybe sse_local)

sse_local=0.0D0  ! 20131209
do i_time_step=1,n_time_steps

  fvec(i_time_step)=0.0D0

  !if( L_GP_print .and. GP_para_flag .and. myid == 1 .and. i_time_step == n_time_steps )then
  !    write(GP_print_unit,'(/A,1x,I6, 1x,I10)') &
  !      'fcn: myid, i_time_step', myid, i_time_step
  !endif ! L_GP_print
  !if( L_ga_print )then
  !    write(GA_print_unit,'(/A,1x,I6, 1x,I10)')'fcn: myid, i_time_step ', myid, i_time_step
  !endif ! L_ga_print

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
      !if( L_ga_print )then
      !    write(GA_print_unit,'(A,2(1x,I6), 3(1x,E15.7))') &
      !          'fcn: myid, i_eqn, RK_soln, data_array, var ', &
      !                myid, i_CODE_equation,                   &
      !                Runge_Kutta_Solution(i_time_step,i_CODE_equation), &
      !                Data_Array(i_time_step,i_CODE_equation), &
      !                data_variance(i_CODE_equation)

      !    write(GA_print_unit,'(A,2(1x,I6), 1x,E15.7)') &
      !          'fcn: myid, i_eqn, data_variance ', &
      !                myid, i_CODE_equation, data_variance(i_CODE_equation)
      !endif ! L_ga_print

      !fvec(i_time_step) = fvec(i_time_step)  +                             &
      ! (   Data_Array(i_time_step,i_CODE_equation) -                       &
      !     Runge_Kutta_Solution(i_time_step,i_CODE_equation)   )**2  /     &
      !                                     Data_Variance(i_CODE_equation)

      fvec(i_time_step) = fvec(i_time_step)  +                                &
          (   Data_Array(i_time_step,i_CODE_equation) -                       &
              Runge_Kutta_Solution(i_time_step,i_CODE_equation)   )**2  *     &
                                            Data_Variance_inv(i_CODE_equation)

  enddo ! i_CODE_equation


  sse_local = sse_local + fvec(i_time_step)  ! 20131209

  !if( L_GP_print .and. GP_para_flag .and. myid == 1 .and. i_time_step == n_time_steps )then
  !   write(GP_print_unit,'(A,1x,I6, 1x,I6, 1x, E15.7)')&
  !          'fcn: myid, i_time_step, fvec', &
  !                myid, i_time_step, fvec(i_time_step)
  !endif ! L_GP_print
  !if( L_ga_print )then
  !    write(GA_print_unit,'(A,1x,I6, 1x,I6, 1x, E15.7)')&
  !          'fcn: myid, i_time_step, fvec ', &
  !                myid, i_time_step, fvec(i_time_step)
  !endif ! L_ga_print

enddo ! i_time_step

!if( L_GP_print .and. GP_para_flag .and. myid == 3  )then
!    write(GP_print_unit,'(A,1x,I6,2x,E24.16)') 'fcn: myid, sse_local = ',myid, sse_local
!endif ! L_GP_print
!if( L_ga_print .and. myid == 1 )then
!    write(GA_print_unit,'(A,1x,I6,2x,E24.16)') 'fcn: myid, sse_local = ',myid, sse_local
!endif ! L_ga_print


return


end subroutine fcn
