subroutine setup_run_fcn( i_GA_indiv,  child_parameters, individual_quality )

! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module


use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module
use GP_data_module


implicit none


integer, intent(in)  ::  i_GA_indiv

integer ::  iflag


! lmdif arrays and variables

real (kind=8) :: x_LMDIF(n_maximum_number_parameters)
real (kind=8) :: fvec(n_time_steps)


real(kind=8), parameter :: tol = 1.0d-30
real(kind=8), parameter :: epsfcn = 1.0d-6    ! original
real(kind=8), parameter :: factor=1.0D+0
real(kind=8), parameter :: zero = 0.0d0

integer (kind=4) ::   info


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=4) :: individual_quality(n_GA_individuals)


real(kind=8) :: child_parameters(n_GA_individuals,n_maximum_number_parameters)

external :: fcn

!--------------------------------------------------------------------------------------------

!write(GA_print_unit,'(/A,2(1x,I6)/)') 'setrf: myid, i_GA_indiv ', myid, i_GA_indiv

x_LMDIF(1:n_maximum_number_parameters) = 0.0D0

do  i_parameter=1,n_parameters

    X_LMDIF(i_parameter) = child_parameters(i_GA_indiv,i_parameter)

    !if( myid == 1 )then
    !    write(GA_print_unit,'(A,3(1x,I6),1x,E20.10)') &
    !          'setrf:1 myid, i_GA_indiv,i_parameter, child_parameters ', &
    !                   myid, i_GA_indiv,i_parameter, &
    !                   child_parameters(i_GA_indiv,i_parameter)
    !    write(GA_print_unit,'(A,3(1x,I6),1x,E20.10)') &
    !          'setrf:1 myid, i_GA_indiv, i_parameter,  X_LMDIF', &
    !                   myid, i_GA_indiv, i_parameter,  X_LMDIF(i_parameter)
    !endif ! myid == 1 

enddo ! i_parameter


!----------------------------------------------------------------------------------------

! call fcn


! re-do is turned off
! if you get a bad result from fcn,
! set the child_parameters to a new set of random numbers, and then try again.
! Do this a maximum of n_redo_max times.
! If you get n_redo_max bad results,
! set the child_parameters to -1.0 and go to the next individual


!if( myid == 1 )then
!    write(GA_print_unit,'(/A,4(1x,I10))') &
!          'setrf: call fcn, myid, i_GA_indiv, n_time_steps, n_parameters     ', &
!                            myid, i_GA_indiv, n_time_steps, n_parameters
!endif ! myid == 1 


iflag = 1

call fcn( n_time_steps, n_parameters, x_LMDIF, fvec, iflag )
info = iflag


!if( myid == 1 )then
!    write(GA_print_unit,'(A,5(1x,I10)/)') &
!              'setrf: aft call fcn myid, i_GA_indiv, n_time_steps, n_parameters, info ', &
!                                   myid, i_GA_indiv, n_time_steps, n_parameters, info
!endif ! myid == 1 


!----------------------------------------------------------------------------------------

! if info < 0 , delete this individual

if( info < 0 ) then

    individual_quality( i_GA_indiv ) = -1
    individual_SSE(i_GA_indiv) =  1.0D+12

    write(GA_print_unit,'(A, 3(1x, I6),  1x,E15.7/)') &
          'setrf:3 myid, i_GA_indiv, quality, SSE ', &
                   myid, i_GA_indiv, &
                   individual_quality(i_GA_indiv), &
                   individual_SSE(i_GA_indiv)
    return

endif ! info < 0


!!! if info < 0, the function gave a bad result for this individual
!!! so retry for a maximum of 10 times to see if a good result can be obtained
!!if( info < 0 ) then
!!
!!    n_retry = 0
!!
!!    retry_loop:&
!!    do
!!
!!        n_retry = n_retry + 1
!!
!!        if( n_retry > 10 )then
!!            ! if info < 0, the function gave a bad result for this individual
!!            ! so mark this individual to skip it in later calculations
!!            individual_quality( i_GA_indiv ) = -1
!!            exit retry_loop
!!        endif
!!
!!        do i_parameter=1,n_parameters
!!
!!           call random_real(cff)
!!
!!            if( cff < 0.0d0 ) then                  ! debug only
!!                cff = max( -1.0e6 , cff )           ! debug only
!!            else                                    ! debug only
!!                cff = min(  1.0e6 , cff )           ! debug only
!!            endif                                   ! debug only
!!
!!
!!            child_parameters(i_GA_indiv,i_parameter)=cff
!!            X_LMDIF(i_parameter)=dble(child_parameters(i_GA_indiv,i_parameter))
!!
!!            !write(GA_print_unit,'(A,3(1x,I6),1x,E20.10)') &
!!            !'setrf:2 myid, i_GA_indiv, i_parameter, child_parameters(i_GA_indiv,i_parameter) ', &
!!            !         myid, i_GA_indiv, i_parameter, child_parameters(i_GA_indiv,i_parameter)
!!            !write(GA_print_unit,'(A,2(1x,I6),1x,E20.10)') &
!!            !      'setrf:2 myid, i_parameter,  X_LMDIF(i_parameter) ', &
!!            !               myid, i_parameter,  X_LMDIF(i_parameter)
!!
!!        enddo ! i_parameter
!!
!!
!!         iflag = 1
!!         call fcn( n_time_steps, n_parameters, x_LMDIF, fvec, iflag)
!!
!!         info = iflag
!!
!!        write(GA_print_unit,'(A,1x,3(1x,I6))')&
!!        'setrf: aft call fcn  myid, n_retry, info  ', myid, n_retry, info
!!
!!        if( info >= 0 ) exit retry_loop
!!
!!    enddo retry_loop    ! n_retry
!!
!!endif ! info < 0


if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------

!write(GA_print_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!      'setrf:3 myid, i_GA_indiv, X_LMDIF', &
!               myid, i_GA_indiv, X_LMDIF(1:n_parameters)


do  i_parameter=1,n_parameters
    child_parameters(i_GA_indiv,i_parameter) = &
                            dabs( x_LMDIF(i_parameter) )
    !write(GA_print_unit,'(A,3(1x,I6),1x,E20.10)') &
    !      'setrf:3 aft RK myid, i_GA_indiv, i_parameter,  X_LMDIF', &
    !                      myid, i_GA_indiv, i_parameter,  X_LMDIF(i_parameter)
enddo ! i_parameter

!write(GA_print_unit,'(/A/ 2(1x, I6), 12( 1x,E24.16))') &
!      'setrf:4 myid, i_GA_indiv, child_parameters(i_GA_indiv,:)', &
!               myid, i_GA_indiv, child_parameters(i_GA_indiv,1:n_parameters)

!-----------------------------------------------------------------------------------

!  calculate the individual SSE values by summing fvec over all time steps
!  fvec(i) = ( fcn(i) - truth(i) )**2
!  so SSE is calculated by summing fvec, not fvec**2

!write(GA_print_unit,'(/A/)')'setrf: calculate the individual SSE values '


individual_SSE(i_GA_indiv)=0.0D+0

if( individual_quality( i_GA_indiv ) > 0 ) then

    !write(20,'(A,1x,I6)') 'setrf: i_GA_indiv ', i_GA_indiv
    !write(GA_print_unit,'(A,1x,I6)') 'setrf: i_GA_indiv ', i_GA_indiv

    do i_time_step=1,n_time_steps

       if( isnan(fvec(i_time_step)) )    fvec(i_time_step) = 0.0d0
       if( abs(fvec(i_time_step)) >  1.0d20 ) fvec(i_time_step) =  1.0d20

       !write(20,'(A,1x,I6,1x,E24.16)' ) &
       !      'setrf: i_time_step, fvec(i_time_step) ', &
       !              i_time_step, fvec(i_time_step)
       !write(GA_print_unit,'(A,1x,I6,1x,E24.16)' ) &
       !      'setrf: i_time_step, fvec(i_time_step) ', &
       !              i_time_step, fvec(i_time_step)

       individual_SSE(i_GA_indiv) = individual_SSE(i_GA_indiv) + &
                                    fvec(i_time_step)

    enddo ! i_time_step

endif !  individual_quality( i_GA_indiv ) > 0

!write(GA_print_unit,'(A,3(1x,I6), 1x, E24.16)') &
!      'setrf: myid, i_GA_indiv, individual_quality, individual_SSE', &
!              myid, i_GA_indiv, &
!              individual_quality( i_GA_indiv ), &
!              individual_SSE(i_GA_indiv)



return

end subroutine setup_run_fcn
