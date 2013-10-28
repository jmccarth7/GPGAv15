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



integer(kind=4) ::   info
integer(kind=4) :: i_time_step
integer(kind=4) :: i_parameter


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=4) :: individual_quality(n_GA_individuals)


real(kind=8) :: child_parameters(n_maximum_number_parameters, n_GA_individuals)

external :: fcn


!--------------------------------------------------------------------------------------------

!if( L_ga_print )then                                                                                
!    write(GA_print_unit,'(/A,2(1x,I6)/)') &
!          'setrf: myid, i_GA_indiv ', myid, i_GA_indiv
!endif ! L_ga_print  

x_LMDIF(1:n_maximum_number_parameters) = 0.0D0

do  i_parameter=1,n_parameters

    X_LMDIF(i_parameter) = child_parameters(i_parameter,i_GA_indiv)

    !if( myid == 1 )then
    !    if( L_ga_print )then   
    !        write(GA_print_unit,'(A,3(1x,I6),1x,E20.10)') &
    !              'setrf:1 myid, i_GA_indiv,i_parameter, child_parameters ', &
    !                       myid, i_GA_indiv,i_parameter, &
    !                       child_parameters(i_parameter,i_GA_indiv)
    !        write(GA_print_unit,'(A,3(1x,I6),1x,E20.10)') &
    !              'setrf:1 myid, i_GA_indiv, i_parameter,  X_LMDIF', &
    !                       myid, i_GA_indiv, i_parameter,  X_LMDIF(i_parameter)
    !    endif ! L_ga_print  
    !endif ! myid == 1

enddo ! i_parameter


!----------------------------------------------------------------------------------------

! call fcn


!if( myid == 1 )then
!    if( L_ga_print )then                                                                                
!        write(GA_print_unit,'(/A,4(1x,I10))') &
!              'setrf: call fcn, myid, i_GA_indiv, n_time_steps, n_parameters', &
!                                myid, i_GA_indiv, n_time_steps, n_parameters
!    endif ! L_ga_print  
!endif ! myid == 1


iflag = 1



call fcn( n_time_steps, n_parameters, x_LMDIF, fvec, iflag )



info = iflag


!if( myid == 1 )then
!    if( L_ga_print )then                                                                                
!        write(GA_print_unit,'(A,5(1x,I10)/)') &
!         'setrf: aft call fcn myid, i_GA_indiv, n_time_steps, n_parameters, info ', &
!                              myid, i_GA_indiv, n_time_steps, n_parameters, info
!    endif ! L_ga_print  
!endif ! myid == 1


!----------------------------------------------------------------------------------------

! if info < 0 , delete this individual

if( info < 0 ) then

    individual_quality( i_GA_indiv ) = -1
    individual_SSE(i_GA_indiv) =  1.0D+12

    !if( L_ga_print )then 
    !    write(GA_print_unit,'(A, 3(1x, I6),  1x,E15.7/)') &
    !          'setrf:3 myid, i_GA_indiv, quality, SSE ', &
    !                   myid, i_GA_indiv, &
    !                   individual_quality(i_GA_indiv), &
    !                   individual_SSE(i_GA_indiv)
    !endif ! L_ga_print  
    return

endif ! info < 0



if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------

!if( L_ga_print )then                                                                                
!    write(GA_print_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!          'setrf:3 myid, i_GA_indiv, X_LMDIF', &
!                   myid, i_GA_indiv, X_LMDIF(1:n_parameters)
!endif ! L_ga_print  


do  i_parameter=1,n_parameters

    child_parameters(i_parameter,i_GA_indiv) = &
                            dabs( x_LMDIF(i_parameter) )

    !if( L_ga_print )then  
    !    write(GA_print_unit,'(A,3(1x,I6),1x,E20.10)') &
    !          'setrf:3 aft RK myid, i_GA_indiv, i_parameter,  X_LMDIF', &
    !                          myid, i_GA_indiv, i_parameter,  X_LMDIF(i_parameter)
    !endif ! L_ga_print  

enddo ! i_parameter

!if( L_ga_print )then                                                                                
!    write(GA_print_unit,'(/A/ 2(1x, I6), 12( 1x,E24.16))') &
!          'setrf:4 myid, i_GA_indiv, child_parameters(:,i_GA_indiv)', &
!                   myid, i_GA_indiv, child_parameters(1:n_parameters,i_GA_indiv)
!endif ! L_ga_print  

!-----------------------------------------------------------------------------------

!  calculate the individual SSE values by summing fvec over all time steps
!  fvec(i) = ( fcn(i) - truth(i) )**2
!  so SSE is calculated by summing fvec, not fvec**2

!if( L_ga_print )then                                                                                
!    write(GA_print_unit,'(/A/)')'setrf: calculate the individual SSE values '
!endif ! L_ga_print  


individual_SSE(i_GA_indiv)=0.0D+0

if( individual_quality( i_GA_indiv ) > 0 ) then

    !if( L_ga_print )then   
    !    write(GA_print_unit,'(A,1x,I6)') 'setrf: i_GA_indiv ', i_GA_indiv
    !endif ! L_ga_print  

    do i_time_step=1,n_time_steps

!old       if( isnan(fvec(i_time_step)) )    fvec(i_time_step) = 0.0d0
!old       if( abs(fvec(i_time_step)) >  1.0d20 ) fvec(i_time_step) =  1.0d20

       if( isnan(fvec(i_time_step)) )    fvec(i_time_step) = 0.0d0
       if( abs(fvec(i_time_step)) >  1.0d20 ) fvec(i_time_step) =  1.0d20

!new        if( isnan(fvec(i_time_step))  .or.   &
!new            abs(fvec(i_time_step)) >  1.0d20   ) fvec(i_time_step) =  1.0d20

       !if( L_ga_print )then
       !    write(GA_print_unit,'(A,1x,I6,1x,E24.16)' ) &
       !          'setrf: i_time_step, fvec(i_time_step) ', &
       !                  i_time_step, fvec(i_time_step)
       !endif ! L_ga_print  

       individual_SSE(i_GA_indiv) = individual_SSE(i_GA_indiv) + &
                                    fvec(i_time_step)

    enddo ! i_time_step

endif !  individual_quality( i_GA_indiv ) > 0

!if( L_ga_print )then                                                                                
!    write(GA_print_unit,'(A,3(1x,I6), 1x, E24.16)') &
!          'setrf: myid, i_GA_indiv, individual_quality, individual_SSE', &
!                  myid, i_GA_indiv, &
!                  individual_quality( i_GA_indiv ), &
!                  individual_SSE(i_GA_indiv)
!endif ! L_ga_print  



return

end subroutine setup_run_fcn
