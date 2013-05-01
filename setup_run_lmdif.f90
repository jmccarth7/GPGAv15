subroutine setup_run_lmdif( i_GA_indiv,  child_parameters, individual_quality )

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


! lmdif arrays and variables

real (kind=8) :: x_LMDIF(n_maximum_number_parameters)
real (kind=8) :: fvec(n_time_steps)
real (kind=8) :: ftol,xtol,gtol


real(kind=8), parameter :: tol = 1.0d-30
real(kind=8), parameter :: epsfcn = 1.0d-6    ! original
real(kind=8), parameter :: factor=1.0D+0
real(kind=8), parameter :: zero = 0.0d0

real (kind=8) :: diag(n_parameters)
real (kind=8) :: fjac(n_time_steps,n_parameters)
real (kind=8) :: qtf(n_parameters)
integer (kind=4) :: maxfev,ldfjac,mode,nprint,info,nfev
integer (kind=4) :: ipvt(n_parameters)


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=4) :: individual_quality(n_GA_individuals)


real(kind=8) :: child_parameters(n_GA_individuals,n_maximum_number_parameters)

external :: fcn


integer(kind=4) :: i_GP_individual                                                                           
integer(kind=4) :: i_GP_Best_Parent                                                                          
integer(kind=4) :: i_GP_Generation                                                                           
integer(kind=4) :: j_GP_Individual                                                                           
integer(kind=4) :: i_Tree                                                                                    
integer(kind=4) :: i_Node            

!--------------------------------------------------------------------------------------------


do  i_parameter=1,n_parameters

    X_LMDIF(i_parameter) = child_parameters(i_GA_indiv,i_parameter)

    !write(GA_print_unit,'(A,3(1x,I6),1x,E20.10)') &
    !      'setrlm:1 myid, i_GA_indiv,i_parameter, child_parameters ', &
    !                myid, i_GA_indiv,i_parameter, &
    !                child_parameters(i_GA_indiv,i_parameter)
    !write(GA_print_unit,'(A,2(1x,I6),1x,E20.10)') &
    !      'setrlm:1 myid, i_parameter,  X_LMDIF', &
    !                myid, i_parameter,  X_LMDIF(i_parameter)

enddo ! i_parameter

!write(GA_print_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!      'setrlm:1 myid, i_GA_indiv, X_LMDIF', &
!                myid, i_GA_indiv, X_LMDIF(1:n_parameters)


! for each of these first individuals, optimize the variables using lmdif.f

info = 0

! maximum iterations in lmdif for function evaluation

!off      maxfev=100*(n_time_steps+1)*100

maxfev= 1000 ! 50 ! 10 ! 10000

ftol=1.0D-10
xtol=1.0D-10

gtol=zero

mode=1
info=1

! nprint < 0  means no printout
nprint= 1  ! set back to zero after diag

ldfjac=n_time_steps

!write(GA_print_unit,*) 'setrlm: i_GA_indiv ', i_GA_indiv

!----------------------------------------------------------------------------------------

if( Lprint_lmdif )then
    write(GA_print_unit,'(/A,4(1x,I6))') &
          'setrlm: call lmdif, myid, n_time_steps, n_parameters, i_GA_indiv ', &
                               myid, n_time_steps, n_parameters, i_GA_indiv
    write(GA_print_unit,'(/A)') 'setrlm: lmdif parameters '

    write(GA_print_unit,'(A,3(1x,E15.7))') 'setrlm: ftol, xtol, gtol     ', ftol, xtol, gtol
    write(GA_print_unit,'(A,3(1x,I10))')   'setrlm: mode, nprint, ldfjac ', mode, nprint, ldfjac
    write(GA_print_unit,'(A,3(1x,E15.7))') 'setrlm: tol,epsfcn, factor   ', tol, epsfcn,factor

    write(GA_print_unit,'(A,1x,I10)')   'setrlm: maxfev ', maxfev
    write(GA_print_unit,'(A,1x,I10)')   'setrlm: info   ', info
endif ! Lprint_lmdif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


L_bad_result = .false.

call lmdif( fcn, n_time_steps, n_parameters, x_LMDIF, fvec, &
            ftol, xtol, gtol, maxfev, epsfcn, &
            diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )

if( Lprint_lmdif )then

    write(GA_print_unit,'(A,4(1x,I10)/)') &
          'setrlm: aft call lmdif, myid, n_time_steps, n_parameters, info ', &
                                   myid, n_time_steps, n_parameters, info

    if( info >= 0 ) then

        write(GA_print_unit,'(A,1x,I10/)') 'setrlm: info flag =  ', info

        write(GA_print_unit,'(A/)') &
        '#######################################################################################'
        write(GA_print_unit,'(A)') 'INFO, error flag.  '

        write(GA_print_unit,'(/A)') &
        'If the user has terminated execution, INFO is set to the (negative) value of IFLAG.'
        write(GA_print_unit,'(A)') 'See the description  of FCN.'

        write(GA_print_unit,'(/A/)') 'Otherwise, INFO is set as follows:'

        write(GA_print_unit,'(A)')  '0, improper input parameters.'
        write(GA_print_unit,'(A)')  &
        '1, both actual and predicted relative reductions in the sum of squares are at most FTOL.'
        write(GA_print_unit,'(A)')  &
        '2, relative error between two consecutive iterates is at most XTOL.'
        write(GA_print_unit,'(A)')  '3, conditions for INFO = 1 and INFO = 2 both hold.'
        write(GA_print_unit,'(A)')  '4, the cosine of the angle between FVEC and &
              &any column of the Jacobian is at most GTOL in absolute value.'
        write(GA_print_unit,'(A)')  '5, number of calls to FCN has reached or exceeded MAXFEV.'
        write(GA_print_unit,'(A)')  &
        '6, FTOL is too small.  No further reduction in the sum of squares is possible.'
        write(GA_print_unit,'(A)')  &
        '7, XTOL is too small.  No further improvement in the approximate solution X is possible.'
        write(GA_print_unit,'(A)') '8, GTOL is too small.  FVEC is orthogonal &
              &to the columns of the Jacobian to machine precision.'
        write(GA_print_unit,'(/A/)') &
        '#######################################################################################'

    endif ! info > 0

    Lprint_lmdif = .FALSE.
endif ! Lprint_lmdif
!----------------------------------------------------------------------------------------


! if info < 0 , delete this individual

if( info < 0 ) then

    individual_quality( i_GA_indiv ) = -1
    individual_SSE(i_GA_indiv) =  1.0D+12
    write(GA_print_unit,'(/A/ 3(1x, I6),  1x,E12.5)') &
          'setrlm:3 myid, i_GA_indiv, individual_quality(i_GA_indiv), &
                                      &individual_SSE(i_GA_indiv) ', &
                    myid, i_GA_indiv, individual_quality(i_GA_indiv), &
                                       individual_SSE(i_GA_indiv)
    return

endif ! info < 0




if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------

!write(GA_print_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!      'setrlm:3 myid, i_GA_indiv, X_LMDIF', &
!                myid, i_GA_indiv, X_LMDIF(1:n_parameters)


do  i_parameter=1,n_parameters
    child_parameters(i_GA_indiv,i_parameter) = &
                           dabs( x_LMDIF(i_parameter) )
enddo ! i_parameter

!write(GA_print_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!      'setrlm:4 myid, i_GA_indiv, child_parameters(i_GA_indiv,:)', &
!                myid, i_GA_indiv, child_parameters(i_GA_indiv,1:n_parameters)


!-----------------------------------------------------------------------------------


!  calculate the individual SSE values by summing fvec over all time steps
!  fvec(i) = ( fcn(i) - truth(i) )**2
!  so SSE is calculated by summing fvec, not fvec**2



!write(GA_print_unit,'(/A/)')'setrlm: calculate the individual SSE values '


if( individual_quality( i_GA_indiv ) > 0 ) then

    !write(10,*) 'setrlm: i_GA_indiv ', i_GA_indiv

    individual_SSE(i_GA_indiv)=0.0D+0

    do i_time_step=1,n_time_steps

       if( isnan(fvec(i_time_step)) ) fvec(i_time_step) = 0.0d0
       if( abs(fvec(i_time_step)) >  1.0d20 ) fvec(i_time_step) =  1.0d20

       !write(10, *) 'setrlm: i_time_step, fvec(i_time_step) ', &
       !                      i_time_step, fvec(i_time_step)

       individual_SSE(i_GA_indiv) = individual_SSE(i_GA_indiv) + fvec(i_time_step)

    enddo ! i_time_step

endif !  individual_quality( i_GA_indiv ) > 0

!write(GA_print_unit,'(A,3(1x,I6), 1x, E15.7)') &
!      'setrlm: myid, i_GA_indiv, individual_quality, individual_SSE', &
!               myid, i_GA_indiv, &
!               individual_quality( i_GA_indiv ), &
!               individual_SSE(i_GA_indiv)



return

end subroutine setup_run_lmdif
