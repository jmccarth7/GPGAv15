subroutine setup_run_para_lmdif( i_G_indiv,  child_parameters, &
                                 individual_quality, &
                                 n_indiv, my_indiv_SSE, &
                                 n_parms, &
                                 L_myprint, myprint_unit  )

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


integer, intent(in)  ::  i_G_indiv
integer, intent(in)  ::  n_indiv
integer, intent(in)  ::  n_parms

real(kind=8),dimension(n_indiv)  ::  my_indiv_SSE

logical, intent(in)  ::  L_myprint
integer, intent(in)  ::  myprint_unit

! lmdif arrays and variables

real(kind=8) :: x_LMDIF(n_maximum_number_parameters)
real(kind=8) :: fvec(n_time_steps)
real(kind=8) :: ftol,xtol,gtol


real(kind=8), parameter :: tol = 1.0d-30
real(kind=8), parameter :: epsfcn = 1.0d-9   ! 1.0d-6    ! original
real(kind=8), parameter :: factor=1.0D+0
real(kind=8), parameter :: zero = 0.0d0

real(kind=8) :: diag(n_maximum_number_parameters)
real(kind=8) :: fjac(n_time_steps,n_maximum_number_parameters)
real(kind=8) :: qtf(n_maximum_number_parameters)
integer(kind=4) :: maxfev,ldfjac,mode,nprint,info,nfev
integer(kind=4) :: ipvt(n_maximum_number_parameters)


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=4) :: individual_quality(n_indiv)

integer(kind=4) :: i_time_step
integer(kind=4) :: i_parameter

real(kind=8) :: child_parameters( n_maximum_number_parameters, n_indiv )

external :: fcn

real(kind=8) :: t1
real(kind=8) :: t2
!real(kind=8) :: delta_wt

!--------------------------------------------------------------------------------------------

!!!!!!!!n_parameters = GP_n_parms( i_G_indiv )

!write(myprint_unit,'(//A,3(1x,I6),1x,E20.10)') &
!      'strplm:1 at entry myid, myprint_unit, n_parms', &
!                         myid, myprint_unit, n_parms
!write(myprint_unit,'(/A,1x,I3,2(1x,I10))') &
!      'strplm: at entry myid, i_G_indiv, individual_quality(i_G_indiv)', &
!                        myid, i_G_indiv, individual_quality(i_G_indiv)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

GP_Individual_Node_Type(1:n_Nodes,1:n_Trees) = &
               GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,i_G_indiv)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

do  i_parameter=1,n_parms

    X_LMDIF(i_parameter) = child_parameters(i_parameter, i_G_indiv)

    !if( L_myprint )then
        !write(myprint_unit,'(A,3(1x,I6),1x,E20.10)') &
        !  'strplm:1 myid, i_G_indiv,i_parameter, child_parameters', &
        !            myid, i_G_indiv,i_parameter, &
        !            child_parameters(i_parameter, i_G_indiv)
        !write(myprint_unit,'(A,2(1x,I6),1x,E20.10)') &
        !  'strplm:1 myid, i_parameter,  X_LMDIF', &
        !            myid, i_parameter,  X_LMDIF(i_parameter)
    !endif ! L_myprint

enddo ! i_parameter

!if( L_myprint )then
!    write(myprint_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!          'strplm:1 myid, i_G_indiv, X_LMDIF', &
!                    myid, i_G_indiv, X_LMDIF(1:n_parms)
!endif ! L_myprint


! for each of these individuals, optimize the variables using lmdif.f

info = 0

! maximum iterations in lmdif for function evaluation

!off      maxfev=100*(n_time_steps+1)*100

maxfev= 2000 ! 50 ! 10 ! 10000

ftol=1.0D-10
xtol=1.0D-10

gtol=zero

mode=1
info=0 ! 1

! nprint < 0  means no printout
nprint= 1  ! set back to zero after diag

ldfjac = n_time_steps

!if( L_myprint )then
!    write(myprint_unit,'(A,1x,I10)') &
!         'strplm: i_G_indiv', i_G_indiv
!endif ! L_myprint

!----------------------------------------------------------------------------------------

!if( Lprint_lmdif )then
!    if( L_myprint )then
!        !write(myprint_unit,'(/A,4(1x,I6))') &
!        ! 'strplm: call lmdif, myid, n_time_steps, n_parms, i_G_indiv', &
!        !                      myid, n_time_steps, n_parms, i_G_indiv
!        !write(myprint_unit,'(/A)') 'strplm: lmdif parameters '
!
!        !write(myprint_unit,'(A,3(1x,I10))')   'strplm: mode, nprint, ldfjac', &
!        !                                               mode, nprint, ldfjac
!        !write(myprint_unit,'(A,3(1x,E15.7))') 'strplm: ftol, xtol, gtol    ', &
!        !                                               ftol, xtol, gtol
!        !write(myprint_unit,'(A,3(1x,E15.7))') 'strplm: tol,epsfcn, factor  ', &
!        !                                               tol, epsfcn,factor
!        !write(myprint_unit,'(A,1x,I10)')   'strplm: maxfev', maxfev
!        !write(myprint_unit,'(A,1x,I10)')   'strplm: info  ', info
!    endif ! L_myprint
!endif ! Lprint_lmdif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


L_bad_result = .false.

t1 = MPI_Wtime()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!write(myprint_unit,'(/A,1x,I3/)') 'strplm: RUN LMDIF myid =', myid


call lmdif( fcn, n_time_steps, n_parms, x_LMDIF, fvec, &
            ftol, xtol, gtol, maxfev, epsfcn, &
            diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

t2 = MPI_Wtime()
!delta_wt = MPI_Wtick()

sum_lmdif = sum_lmdif + ( t2 - t1 )

!write(myprint_unit,'(/A,1x,E15.7)') 'strplm: time spent in lmdif = ', t2 - t1
!write(myprint_unit,'(A,1x,E15.7)')  'strplm: time increment      = ', delta_wt
!write(myprint_unit,'(A,1x,E15.7/)') 'strplm:          sum_ lmdif = ', sum_lmdif

!write(6,'(A,3(1x,I3),1x,I10/)') &
!      'strplm: aft call lmdif, myid, n_parms, info, n_time_steps', &
!                               myid, n_parms, info, n_time_steps

if( Lprint_lmdif )then

    if( L_myprint )then

        !write(myprint_unit,'(A,3(1x,I3),1x,I10/)') &
        !      'strplm: aft call lmdif, myid, n_parms, info, n_time_steps', &
        !                               myid, n_parms, info, n_time_steps

        !!if( info >= 0 ) then
        !!
        !!    write(myprint_unit,'(A,1x,I10/)') 'strplm: info flag =  ', info
        !!
        !!    write(myprint_unit,'(A/)') &
        !!    '######################################################################################'
        !!    write(myprint_unit,'(A)') 'INFO, error flag.  '
        !!
        !!    write(myprint_unit,'(/A)') &
        !!    'If the user has terminated execution, &
        !!    &INFO is set to the (negative) value of IFLAG.'
        !!    write(myprint_unit,'(A)') 'See the description  of FCN.'
        !!
        !!    write(myprint_unit,'(/A/)') 'Otherwise, INFO is set as follows:'
        !!
        !!    write(myprint_unit,'(A)')  '0, improper input parameters.'
        !!    write(myprint_unit,'(A)')  &
        !!    '1, both actual and predicted relative reductions &
        !!    &in the sum of squares are at most FTOL.'
        !!    write(myprint_unit,'(A)')  &
        !!    '2, relative error between two consecutive iterates is at most XTOL.'
        !!    write(myprint_unit,'(A)')  '3, conditions for INFO = 1 and INFO = 2 both hold.'
        !!    write(myprint_unit,'(A)')  '4, the cosine of the angle between FVEC and &
        !!          &any column of the Jacobian is at most GTOL in absolute value.'
        !!    write(myprint_unit,'(A)')  &
        !!    '5, number of calls to FCN has reached or exceeded MAXFEV.'
        !!    write(myprint_unit,'(A)')  &
        !!    '6, FTOL is too small.  No further reduction in the sum of squares is possible.'
        !!    write(myprint_unit,'(A)')  &
        !!    '7, XTOL is too small. &
        !!    &No further improvement in the approximate solution X is possible.'
        !!    write(myprint_unit,'(A)') '8, GTOL is too small.  FVEC is orthogonal &
        !!          &to the columns of the Jacobian to machine precision.'
        !!    write(myprint_unit,'(/A/)') &
        !!    '######################################################################################'
        !!
        !!endif ! info > 0
    endif ! L_myprint

    Lprint_lmdif = .FALSE.
endif ! Lprint_lmdif
!----------------------------------------------------------------------------------------


! if info < 0 , delete this individual

if( info < 0 ) then

    individual_quality( i_G_indiv ) = -1
    my_indiv_SSE(i_G_indiv) =  1.0D+12

    if( L_myprint )then
        write(myprint_unit,'(/A/ 3(1x, I6),  1x,E12.5)') &
          'strplm:3 myid, i_G_indiv, individual_quality(i_G_indiv), &
                                      &my_indiv_SSE(i_G_indiv)', &
                    myid, i_G_indiv, individual_quality(i_G_indiv), &
                                       my_indiv_SSE(i_G_indiv)
    endif ! L_myprint
    return

endif ! info < 0




if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------

!if( L_myprint )then
!    write(myprint_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!          'strplm:3 myid, i_G_indiv, X_LMDIF', &
!                    myid, i_G_indiv, X_LMDIF(1:n_parms)
!endif ! L_myprint


do  i_parameter=1,n_parms
    child_parameters(i_parameter,i_G_indiv) = &
                           dabs( x_LMDIF(i_parameter) )
enddo ! i_parameter

!if( L_myprint )then
!    write(myprint_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!     'strplm:4 myid, i_G_indiv, child_parameters(:,i_G_indiv)', &
!               myid, i_G_indiv, child_parameters(1:n_parms, i_G_indiv)
!endif ! L_myprint


!-----------------------------------------------------------------------------------


!  calculate the individual SSE values by summing fvec over all time steps
!  fvec(i) = ( fcn(i) - truth(i) )**2
!  so SSE is calculated by summing fvec, not fvec**2



!if( L_myprint )then
!    write(myprint_unit,'(/A/)')'strplm: calculate the individual SSE values '
!endif ! L_myprint


if( individual_quality( i_G_indiv ) > 0 ) then

    !write(10,*) 'strplm: i_G_indiv ', i_G_indiv

    my_indiv_SSE(i_G_indiv)=0.0D+0

    do i_time_step=1,n_time_steps

       if( isnan(fvec(i_time_step)) .or.  &
           abs(fvec(i_time_step)) >  1.0d20 ) fvec(i_time_step) =  1.0d20

       !write(10, *) 'strplm: i_time_step, fvec(i_time_step) ', &
       !                      i_time_step, fvec(i_time_step)

       my_indiv_SSE(i_G_indiv) = my_indiv_SSE(i_G_indiv) + fvec(i_time_step)

    enddo ! i_time_step

endif !  individual_quality( i_G_indiv ) > 0

!if( L_myprint )then
!    write(myprint_unit,'(A,3(1x,I6), 1x, E15.7)') &
!      'strplm: myid, i_G_indiv, individual_quality, my_indiv_SSE', &
!               myid, i_G_indiv, &
!               individual_quality( i_G_indiv ), &
!               my_indiv_SSE(i_G_indiv)
!endif ! L_myprint



return

end subroutine setup_run_para_lmdif
